#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys, copy, os, importlib
from mpi4py import MPI
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
from pysimulators import BeamGaussian
from pysimulators.interfaces.healpy import SceneHealpixCMB
from pysimulators.interfaces.healpy import HealpixConvolutionGaussianOperator
import simulation.pointing.generate_pointing as gen_p
from simulation.beam.beam_kernel import get_beam

class Bolo:

    def __init__(self, settings, pointing_params, beam_params, bolo_name):
            self.settings = settings
            self.pointing_params = pointing_params
            self.beam_params = beam_params
            self.bolo_params = importlib.import_module("simulation.bolo.bolo_params." + bolo_name).bolo_params 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def simulate_timestream(self, comm, segment, sky_map):
        
        #Getting the beam profile and the del_beta
        beam_kernel, del_beta = get_beam(self.beam_params)
        
        #Building the projection matrix P
        nsamples = int(1000.0*self.pointing_params.t_segment/self.pointing_params.t_sampling)*self.settings.oversampling_rate 
        npix = hp.nside2npix(self.settings.nside_in)
        matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
        P = ProjectionOperator(matrix)

        signal = np.zeros(nsamples)
        matrix.data.value[:, 0, 0, 0] = 0.5

        for i in range(del_beta.size):
            v, pol_ang = gen_p.generate_pointing(segment, self.pointing_params, self.bolo_params, del_beta[i])
            hit_pix = hp.vec2pix(self.settings.nside_in, v[...,0], v[...,1], v[...,2])
            matrix.data.index[:, 0] = hit_pix
            matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol_ang)
            matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol_ang)
            P.matrix = matrix
            if i is del_beta.size/2:
                v_central = v[::self.settings.oversampling_rate]
            #Generating the time ordered signal
            signal += np.convolve(P(sky_map.T), beam_kernel.T[i], mode = 'same')

        beam_sum = np.sum(beam_kernel)
        signal/=beam_sum
        signal = signal[::self.settings.oversampling_rate]

        hitmap = self.get_hitmap(v_central)

        if self.settings.write_signal:
            self.write_ts(signal, self.bolo_params.bolo_name, segment)

        return hitmap


    def get_hitmap(self, v):
        nsamples = v.shape[0]
        npix = hp.nside2npix(self.settings.nside_in)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
        matrix.data.value = 1
        P = ProjectionOperator(matrix, shapein=npix, shapeout=nsamples)        
        hit_pix = hp.vec2pix(self.settings.nside_in, v[...,0], v[...,1], v[...,2])
        matrix.data.index = hit_pix[..., None]
        P.matrix = matrix
        hitmap = P.T(np.ones(nsamples, dtype=np.float32))
        return hitmap


    def write_ts(self, signal, bolo_name, segment): 
        out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp, bolo_name)
        ts_file = "ts_" + str(segment+1).zfill(4)
        np.save(os.path.join(out_dir, ts_file), signal)

def get_scanned_map(sky_map, hitmap):
    valid = hitmap>0
    sky_map[...,~valid] = np.nan
    return sky_map

def create_output_dirs(settings):
    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp)
    for bolo_name in settings.bolo_names:
        os.makedirs(os.path.join(out_dir, bolo_name))

def get_sky_map(settings):
    sky_map = np.array(hp.read_map(settings.input_map, field=(0,1,2)))
    if not settings.do_pol:
        sky_map[1:,...] = 0.0
    return sky_map

def run_serial(settings, pointing_params, beam_params):
    num_segments = int(pointing_params.t_flight/pointing_params.t_segment)
    sky_map = get_sky_map(settings) 
    hitmap = np.zeros(sky_map.size)
    bolo_num = 0
    count = 0

    for bolo_name in settings.bolo_names:
        bolo = Bolo(settings, pointing_params, beam_params, bolo_name)
        print "Doing Bolo : ", bolo_name
        for i in range(num_segments):
            print "Segment : ", i
            print "Rank : ", count
            hitmap += bolo.simulate_timestream(None, count, sky_map)
            count += 1
        bolo_num += 1

    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp)
    if settings.write_scanned_map:
        scanned_map = get_scanned_map(sky_map, hitmap)
        hp.write_map(os.path.join(out_dir, "scanned_map.fits"), scanned_map)
    if settings.write_hitmap:
        hp.write_map(os.path.join(out_dir, "hitmap_in.fits"), hitmap)

def run_mpi(settings, pointing_params, beam_params):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    num_segments = int(pointing_params.t_flight/pointing_params.t_segment)
    sky_map = get_sky_map(settings) 
    if rank is 0:
        create_output_dirs(settings)
        hitmap = np.zeros(hp.nside2npix(settings.nside_in), dtype=np.float32)
    else:
        hitmap = None
    count = 0
    for bolo_name in settings.bolo_names:
        for segment in range(num_segments):
            if count%size is rank:
                print "Doing Bolo : ", bolo_name, "Segment : ", segment, "Rank : ", rank, "Count :", count
                bolo = Bolo(settings, pointing_params, beam_params, bolo_name)
                hitmap_local = bolo.simulate_timestream(comm, segment, sky_map)
                comm.Reduce(hitmap_local, hitmap, MPI.SUM, 0) 
            count += 1

    if rank is 0:
        out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp)
        if settings.write_scanned_map:
            scanned_map = get_scanned_map(sky_map, hitmap)
            hp.write_map(os.path.join(out_dir, "scanned_map.fits"), scanned_map)
        if settings.write_hitmap:
            hp.write_map(os.path.join(out_dir, "hitmap_in.fits"), hitmap)

if __name__=="__main__":
    from custom_settings import settings, pointing_params, beam_params
    run_mpi(settings, pointing_params, beam_params)
