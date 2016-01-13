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

    hit_pix = hp.vec2pix(settings.nside_out, v[...,0], v[...,1], v[...,2])
    nsamples = hit_pix.size
    npix = hp.nside2npix(settings.nside_out)
   
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


def get_signal(settings, bolo_name, segment):
    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.scanning_time_stamp)
    file_suffix = str(segment+1).zfill(4) + '.npy'
    ts = np.load(os.path.join(out_dir, bolo_name, 'ts_'+file_suffix))
    v = np.load(os.path.join(out_dir, bolo_name, 'vec_'+file_suffix))
    return ts, v

def write_map(settings, sky_map, hitmap):
    out_dir = os.path.join(settings.global_output_dir, "reconstructing", settings.time_stamp)
    os.makedirs(out_dir)
    hp.write_map(os.path.join(out_dir, "reconstructed_map.fits"), sky_map)
    hp.write_map(os.path.join(out_dir, "hitmap_out.fits"), hitmap)

def run_mpi(settings):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size() 
    
    num_segments = int(settings.t_flight/settings.t_segment)
    count = 0

    for bolo_name in settings.bolo_names:
        for segment in range(num_segments): 
            signal, v = get_signal(settings, bolo_name, segment)
            if count%size is rank:
                print "Doing Bolo : ", bolo_name, "Segment : ", segment, "Rank : ", rank, "Count : ", count 
                sky_map, hitmap = make_map_from_signal(signal, v, bolo_name, segment)
            count+=1

    if rank is 0:
        write_map(settings, sky_map, hitmap)
        
if __name__=="__main__":
    from custom_settings import settings
    run_mpi(settings)
