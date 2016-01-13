#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpi4py import MPI
from pyoperators import MPIDistributionIdentityOperator
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
import os, sys

def make_map_from_signal(rank, signal, v):

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

def get_signal(rank, settings):
    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.scanning_time_stamp)
    signal_name = "ts_" + str(rank+1).zfill(4) + ".npy"     
    v_name = "vec_" + str(rank+1).zfill(4) + ".npy"
    signal = np.load(os.path.join(out_dir, "bolo_ts", signal_name))
    v = np.load(os.path.join(out_dir, "pointing", v_name))
    return signal, v

def write_map(settings, sky_map, hitmap):
    out_dir = os.path.join(settings.global_output_dir, "reconstructing", settings.time_stamp)
    os.makedirs(out_dir)
    hp.write_map(os.path.join(out_dir, "reconstructed_map.fits"), sky_map)
    hp.write_map(os.path.join(out_dir, "hitmap_out.fits"), hitmap)

def run_mpi(settings):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size() 
    signal, v = get_signal(rank, settings)
    sky_map, hitmap = make_map_from_signal(rank, signal, v)
    if rank is 0:
        write_map(settings, sky_map, hitmap)

if __name__=="__main__":
    from custom_settings import settings
    run_mpi(settings)
    
