#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys, copy, os
#from mpi4py import MPI
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix
from pysimulators import BeamGaussian
from pysimulators.interfaces.healpy import SceneHealpixCMB
from pysimulators.interfaces.healpy import HealpixConvolutionGaussianOperator
import simulation.pointing.generate_pointing as gen_p
from simulation.beam.beam_kernel import get_beam

class Bolo:

    def __init__(self, settings=None, pointing_params=None, beam_params=None, bolo_name='0001'):
        if pointing_params is None:
            from custom_settings import pointing_params 
            self.pointing_params = pointing_params  
        else:
            self.pointing_params = pointing_params

        if beam_params is None:
            from custom_settings import beam_params 
            self.beam_params = beam_params  
        else:
            self.beam_params = beam_params

        if settings is None:
            from custom_settings import settings 
            self.settings = settings
        else:
            self.settings = settings


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def simulate_timestream(self, rank, sky_map=None):
        
        #Getting the beam profile and the del_beta
        beam_kernel, del_beta = get_beam(self.beam_params)
        
        if sky_map is None:
            sky_map = hp.read_map(self.settings.input_map)

        #Building the projection matrix P
        nsamples = int(1000.0*self.pointing_params.t_segment/self.pointing_params.t_sampling)*self.settings.oversampling_rate 
        npix = hp.nside2npix(self.settings.nside_in)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
        matrix.data.value = 1
        P = ProjectionOperator(matrix, shapein=npix, shapeout=nsamples)        

        signal = np.zeros(nsamples)

        for i in range(del_beta.size):
            v = gen_p.generate_pointing(rank, self.pointing_params, np.deg2rad(del_beta[i]/60.0))
            hit_pix = hp.vec2pix(self.settings.nside_in, v[...,0], v[...,1], v[...,2])
            matrix.data.index = hit_pix[..., None]
            P.matrix = matrix
            if i is del_beta.size/2:
                v_central = v[::self.settings.oversampling_rate]
            #Generating the time ordered signal
            signal += np.convolve(P(sky_map), beam_kernel.T[i], mode = 'same')

        beam_sum = np.sum(beam_kernel)
        signal/=beam_sum
        signal = signal[::self.settings.oversampling_rate]

        hitmap = self.get_hitmap(v_central)

        if self.settings.write_signal:
            self.write_ts(rank, signal)

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


    def write_ts(self, rank, signal): 
        out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp, "bolo_ts")
        ts_file = "ts_" + str(rank+1).zfill(4)
        if rank is 0:
            os.makedirs(out_dir)
        np.save(os.path.join(out_dir, ts_file), signal)

def get_scanned_map(sky_map, hitmap):
    valid = hitmap>0
    sky_map[~valid] = np.nan
    return sky_map

def run_serial(settings, pointing_params, beam_params):
    num_segments = int(pointing_params.t_flight/pointing_params.t_segment)
    sky_map = hp.read_map(settings.input_map)
    hitmap = np.zeros(sky_map.size)
    bolo_num = 0
    count = 0

    for bolo_name in settings.bolo_names:
        bolo = Bolo(settings, pointing_params, beam_params, bolo_name)
        print "Doing Bolo : ", bolo_name
        for i in range(num_segments):
            print "Segment : ", i
            print "Rank : ", count
            hitmap += bolo.simulate_timestream(count, sky_map)
            count += 1
        bolo_num += 1

    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp)
    if settings.write_scanned_map:
        scanned_map = get_scanned_map(sky_map, hitmap)
        hp.write_map(os.path.join(out_dir, "scanned_map.fits"), scanned_map)
    if settings.write_hitmap:
        hp.write_map(os.path.join(out_dir, "hitmap_in.fits"), hitmap)

"""
def run_mpi(settings, pointing_params, beam_params):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    num_segments = int(pointing_params.t_flight/pointing_params.t_segment)
    hitmap = np.zeros(sky_map.size)
    count = 1
    
    for bolo_name in settings.bolo_names:
        for i in range(num_segments):
            if count%size is (rank + 1):
                bolo = Bolo(settings, pointing_params, beam_params, bolo_name)
                comm.Reduce(bolo.simulate_timestream(rank), hitmap, MPI.SUM, 0) 
            count += 1

    if rank is 0:
        out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp)
        if settings.write_scanned_map:
            scanned_map = get_scanned_map(sky_map, hitmap)
            hp.write_map(os.path.join(out_dir, "scanned_map.fits"), scanned_map)
        if settings.write_hitmap:
            hp.write_map(os.path.join(out_dir, "hitmap_in.fits"), hitmap)
"""

if __name__=="__main__":
    from custom_settings import settings, pointing_params, beam_params
    run_serial(settings, pointing_params, beam_params)
