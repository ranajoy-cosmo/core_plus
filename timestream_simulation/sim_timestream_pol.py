#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
import copy
import os
import importlib
#import h5py
import time
#from mpi4py import MPI
from pysimulators import ProjectionOperator, BeamGaussian
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
from simulation.beam.beam_kernel_cartesian import get_beam
import simulation.timestream_simulation.sim_pointing as gen_p
from simulation.lib.utilities.time_util import get_time_stamp

class Bolo:

    def __init__(self, bolo_name):
            self.bolo_params = importlib.import_module("simulation.timestream_simulation.bolo_params." + bolo_name).bolo_params

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for a given bolo with any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def simulate_timestream(self, segment, sky_map, segment_group):

        time_segment_start = time.time()
        #Getting the beam profile and the del_beta
        beam_kernel, del_beta = get_beam(beam_params, self.bolo_params)

        #Building the projection matrix P
        nsamples = int(1000.0*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate
        npix = hp.nside2npix(scan_params.nside)
        matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
        P = ProjectionOperator(matrix)

        signal = np.zeros(nsamples)
        matrix.data.value[:, 0, 0, 0] = 0.5

        t_start = scan_params.t_segment*segment

        for i in range(del_beta.size):
            v, pol_ang = gen_p.generate_pointing(scan_params, self.bolo_params, segment_group, t_start, del_beta[i])
            hit_pix = hp.vec2pix(scan_params.nside, v[...,0], v[...,1], v[...,2])
            matrix.data.index[:, 0] = hit_pix
            matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol_ang)
            matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol_ang)
            P.matrix = matrix
            if i is del_beta.size/2:
                np.save(os.path.join(segment_group, "vector"), v)
                np.save(os.path.join(segment_group, "pol_ang"), pol_ang)
                v_central = v[::scan_params.oversampling_rate]
            #Generating the time ordered signal
            signal += np.convolve(P(sky_map.T), beam_kernel[i], mode = 'same')

        beam_sum = np.sum(beam_kernel)
        signal/=beam_sum
        signal = signal[::scan_params.oversampling_rate]

        hitmap = self.get_hitmap(v_central)

        #segment_group.create_dataset("ts_signal", data=signal)
        np.save(os.path.join(segment_group, "ts_signal"), signal)

        time_segment_end = time.time() 
        print "Time taken for scanning segment :", (time_segment_end - time_segment_start), "seconds"

        return hitmap

    def get_hitmap(self, v):
        nsamples = v.shape[0]
        npix = hp.nside2npix(scan_params.nside)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
        matrix.data.value = 1
        P = ProjectionOperator(matrix, shapein=npix, shapeout=nsamples)
        hit_pix = hp.vec2pix(scan_params.nside, v[...,0], v[...,1], v[...,2])
        matrix.data.index = hit_pix[..., None]
        P.matrix = matrix
        hitmap = P.T(np.ones(nsamples, dtype=np.float32))
        return hitmap


def calculate_params():
    if scan_params.mode is 1:
        scan_params.t_spin = 360.0*60.0*np.sin(scan_params.beta)*scan_params.t_sampling/1000.0/scan_params.theta_co
        scan_params.t_prec = 360.0*60.0*np.sin(scan_params.alpha)*scan_params.t_spin/scan_params.theta_cross

    if scan_params.mode is 2:
        scan_params.theta_cross = 360.0*60.0*np.sin(scan_params.alpha)*scan_params.t_spin/scan_params.t_prec
        scan_params.theta_co = 360*60*np.sin(scan_params.beta)*scan_params.t_sampling/1000.0/scan_params.t_spin

    beam_params.beam_resolution = scan_params.theta_co/scan_params.oversampling_rate

def display_params():
    print "alpha : ", scan_params.alpha, " degrees"
    print "beta : ", scan_params.beta, " degrees"
    print "T flight : ", scan_params.t_flight/60.0/60.0, "hours"
    print "T segment :", scan_params.t_segment/60.0/60.0, "hours"
    print "T precession : ", scan_params.t_prec/60.0/60.0, "hours"
    print "T spin : ", scan_params.t_spin, " seconds"
    print "T sampling : ", scan_params.t_sampling, " milli-seconds"
    print "Scan frequency : ", 1000.0/scan_params.t_sampling, "Hz"
    print "Theta co : ", scan_params.theta_co, " arcmin"
    print "Theta cross : ", scan_params.theta_cross, " arcmin"
    n_steps = int(1000*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate
    print "#Samples per segment : ", n_steps
    print "Estimated use of memory : ", 15*n_steps*8.0/1024/1024, "MB"

def get_scanned_map(sky_map, hitmap):
    valid = hitmap>0
    sky_map[...,~valid] = np.nan
    return sky_map

def get_sky_map():
    sky_map = np.array(hp.read_map(scan_params.input_map, field=(0,1,2)))
    if scan_params.do_only_T:
        sky_map[1:,...] = 0.0
    return sky_map


def run_serial():

    num_segments = int(scan_params.t_flight/scan_params.t_segment)
    sky_map = get_sky_map() 
    hitmap = np.zeros(hp.nside2npix(scan_params.nside))

    time_stamp = get_time_stamp()
    out_dir = os.path.join(scan_params.global_output_dir, "scanning", time_stamp)
    os.makedirs(out_dir)
    
    root_file = h5py.File(os.path.join(out_dir, "data.hdf5"), libver="latest")

    calculate_params()

    if scan_params.display_params:
        display_params()

    for bolo_name in scan_params.bolo_names:
        bolo = Bolo(bolo_name)
        bolo_group = root_file.create_group(bolo_name)
        print "Doing Bolo : ", bolo_name
        for segment in range(num_segments):
            segment_name = str(segment+1).zfill(4)
            segment_group = bolo_group.create_group(segment_name)
            print "Segment : ", segment_name
            hitmap += bolo.simulate_timestream(segment, sky_map, segment_group)

    scanned_map = get_scanned_map(sky_map, hitmap)
    hp.write_map(os.path.join(out_dir, "scanned_map.fits"), scanned_map)
    hp.write_map(os.path.join(out_dir, "hitmap_in.fits"), hitmap)

def make_output_dirs(out_dir, bolo_names, num_segments):
    os.makedirs(out_dir)
    for bolo in bolo_names:
        os.makedirs(os.path.join(out_dir, bolo))
        for segment in range(num_segments):
            segment_name = str(segment+1).zfill(4)
            os.makedirs(os.path.join(out_dir, bolo, segment_name))
            

def run_mpi():
    time_start = time.time()   

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    num_segments = int(scan_params.t_flight/scan_params.t_segment)
    sky_map = get_sky_map()

    time_stamp = [None]
    if rank is 0:
        time_stamp[0] = get_time_stamp()
    time_stamp = comm.bcast(time_stamp, root=0)

    out_dir = os.path.join(scan_params.global_output_dir, "scanning", time_stamp[0])

    if rank is 0:
        make_output_dirs(out_dir, scan_params.bolo_names, num_segments)
        hitmap = np.zeros(hp.nside2npix(scan_params.nside), dtype=np.float32)
    else:
        hitmap = None

    comm.Barrier()

    #root_file = h5py.File(os.path.join(out_dir, "data.hdf5"), driver='mpio', libver="latest", comm=comm)

    calculate_params()

    count = 0
    for bolo_name in scan_params.bolo_names:
        for segment in range(num_segments):
            if count%size is rank:
                #bolo_group = root_file.require_group(bolo_name)
                segment_name = str(segment+1).zfill(4)
                #segment_group = bolo_group.create_group(segment_name)
                segment_group = os.path.join(out_dir, bolo_name, segment_name)
                print "Doing Bolo :", bolo_name, " Segment :", segment, " Rank :", rank, " Count :", count
                bolo = Bolo(bolo_name)
                hitmap_local = bolo.simulate_timestream(segment, sky_map, segment_group)
                comm.Reduce(hitmap_local, hitmap, MPI.SUM, 0)
            count += 1

    if rank is 0:
        #out_dir = os.path.join(scan_params.global_output_dir, "scanning", scan_params.time_stamp)
        #print out_dir
        scanned_map = get_scanned_map(sky_map, hitmap)
        hp.write_map(os.path.join(out_dir, "scanned_map.fits"), scanned_map)
        hp.write_map(os.path.join(out_dir, "hitmap_in.fits"), hitmap)

    time_end = time.time() 
    if rank is 0:
        print "Total time taken :", (time_end- time_start), "seconds"
        print "Scan time stamp :", time_stamp[0]

if __name__=="__main__":
    from custom_params import scan_params, beam_params
    calculate_params()
    display_params()
    sys.exit()
    run_serial()

