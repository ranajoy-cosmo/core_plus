#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpi4py import MPI
from pyoperators import DiagonalOperator, PackOperator, pcg, MPIDistributionIdentityOperator
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
import os

def make_map_from_signal(signal, v, bolo_name, segment):

    hit_pix = hp.vec2pix(map_making_params.nside_out, v[...,0], v[...,1], v[...,2])
    nsamples = hit_pix.size
    npix = hp.nside2npix(map_making_params.nside_out)
   
    matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.value = 1
    matrix.data.index = hit_pix[..., None]

    P = ProjectionOperator(matrix, shapein = npix, shapeout = nsamples)
    P.matrix = matrix
    D = MPIDistributionIdentityOperator()
    H = P*D
    sky_map = (H.T*H).I*H.T*signal

    hitmap = H.T(np.ones(nsamples, dtype=np.float32))

    return sky_map, hitmap


def get_signal(out_dir, bolo_name, segment):
    segment_name = str(segment+1).zfill(4)
    data_dir = os.path.join(out_dir, bolo_name, segment_name)
    ts = np.load(os.path.join(data_dir, "ts_signal.npy"))
    v = np.load(os.path.join(data_dir, "vector.npy"))
    return ts, v

def write_map(sky_map, hitmap):
    out_dir = os.path.join(map_making_params.global_output_dir, "reconstructing", map_making_params.time_stamp)
    os.makedirs(out_dir)
    hp.write_map(os.path.join(out_dir, "reconstructed_map.fits"), sky_map)
    hp.write_map(os.path.join(out_dir, "hitmap_out.fits"), hitmap)

def run_mpi():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size() 
    
    num_segments = int(map_making_params.t_flight/map_making_params.t_segment)
    count = 0

    data_dir = os.path.join(map_making_params.global_output_dir, "scanning", map_making_params.scanning_time_stamp)

    for bolo_name in map_making_params.bolo_names:
        for segment in range(num_segments): 
            signal, v = get_signal(data_dir, bolo_name, segment)
            if count%size is rank:
                print "Doing Bolo : ", bolo_name, "Segment : ", segment, "Rank : ", rank, "Count : ", count 
                sky_map, hitmap = make_map_from_signal(signal, v, bolo_name, segment)
            count+=1

    if rank is 0:
        write_map(sky_map, hitmap)
        
if __name__=="__main__":
    from custom_params import map_making_params
    run_mpi()
