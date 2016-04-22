#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpi4py import MPI
from pyoperators import DiagonalOperator, PackOperator, pcg, MPIDistributionIdentityOperator
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
from memory_profiler import profile
import os
import sys
import shutil
import time

#@profile
def make_map_from_signal(signal, v, pol_ang):

    sys.stdout.flush()
    hit_pix = hp.vec2pix(map_making_params.nside_out, v[...,0], v[...,1], v[...,2])
    nsamples = hit_pix.size
    npix = hp.nside2npix(map_making_params.nside_out)
    
    matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.value = 1
    matrix.data.index[:, 0] = hit_pix
    matrix.data.value[:, 0, 0, 0] = 0.5
    matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol_ang)
    matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol_ang)
    
    P = ProjectionOperator(matrix)#, shapein = npix, shapeout = nsamples)
    D = MPIDistributionIdentityOperator()
    H = P*D
     
    hitmap = H.T(np.ones(nsamples, dtype=np.float32))[:, 0]*2
    mask = hitmap>0

    P = P.restrict_new(mask, inplace=True)
    H = P*D
    pack = PackOperator(mask, broadcast='rightward')

    A = H.T*H
    b = H.T*signal
    M = DiagonalOperator(1/hitmap[mask], broadcast='rightward')

    solution = pcg(A, b, M=M, disp=True, tol=1e-4, maxiter=200)
    x = pack.T*solution['x']
    x[hitmap == 0] = np.nan
    del solution
    del mask
    #sky_map = x.T

    return x.T, hitmap

"""
def get_signal(bolo_name, segment, data_root):
    segment_name = str(segment+1).zfill(4)
    ts = data_root[os.path.join(bolo_name, segment_name, "ts_signal"][:]
    v = data_root[os.path.join(bolo_name, segment_name, "vector"][:]
    pol = data_root[os.path.join(bolo_name, segment_name, "pol_ang"][:]
    return ts, v, pol
"""

def get_signal(out_dir, bolo_name_1, bolo_name_2, segment):
    segment_name = str(segment+1).zfill(4)
    data_dir_1 = os.path.join(out_dir, bolo_name_1, segment_name)
    data_dir_2 = os.path.join(out_dir, bolo_name_2, segment_name)
    ts_1 = np.load(os.path.join(data_dir_1, "ts_signal.npy"))
    ts_2 = np.load(os.path.join(data_dir_2, "ts_signal.npy"))
    v = np.load(os.path.join(data_dir, "vector.npy"))
    pol = np.load(os.path.join(data_dir, "pol_ang.npy"))
    return (ts_1 - ts_2), v, pol

def write_map(sky_map, hitmap):
    out_dir = os.path.join(map_making_params.global_output_dir, "reconstructing", map_making_params.time_stamp)
    param_dir = os.path.join(out_dir, "params")
    os.makedirs(out_dir)
    os.makedirs(param_dir)
    shutil.copy("default_params.py", param_dir)
    shutil.copy("custom_params.py", param_dir)
    shutil.copy("map_maker_pol.py", out_dir)
    hp.write_map(os.path.join(out_dir, "reconstructed_map.fits"), sky_map)
    hp.write_map(os.path.join(out_dir, "hitmap_out.fits"), hitmap)

#@profile
def run_mpi():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size() 
    
    #num_segments = int(map_making_params.t_flight/map_making_params.t_segment)
    count = 0

    data_dir = os.path.join(map_making_params.global_output_dir, "scanning", map_making_params.scanning_time_stamp)
    #data_root = h5py.File(os.path.join(out_dir, "data.hdf5"), "r")
    start = time.time()
    bolo_1 = map_making_params.bolo_names[0]
    bolo_2 = map_making_params.bolo_names[1]

    for segment in range(num_segments): 
        if count%size is rank:
            signal_diff, v, pol_ang = get_diff_signal(data_dir, bolo_1, bolo_2, segment)
            print "Doing Bolos :", bolo_1, "and", bolo_2" Segment :", segment, " Rank :", rank, " Count :", count 
            print signal.size, v.shape, pol_ang.shape
            sky_map, hitmap = make_map_from_signal(signal_diff, v, pol_ang)
        count+=1

    stop = time.time()
    if rank is 0:
        write_map(sky_map, hitmap)
        print "Total time taken :", (stop - start)
        
if __name__=="__main__":
    from custom_params import map_making_params
    run_mpi()
