#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from mpi4py import MPI
from pyoperators import DiagonalOperator, PackOperator, pcg, MPIDistributionIdentityOperator
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
import simulation.bolo.timestream_pol as ts
import os

def make_map_from_signal(signal, v, pol_ang, bolo_name, segment):

    hit_pix = hp.vec2pix(settings.nside_out, v[...,0], v[...,1], v[...,2])
    nsamples = hit_pix.size
    npix = hp.nside2npix(settings.nside_out)
    
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

    H = H.restrict_new(mask, inplace=True)
    pack = PackOperator(mask, broadcast='rightward')

    A = H.T*H
    b = H.T*signal
    M = DiagonalOperator(1/hitmap[mask], broadcast='rightward')
    #M = DiagonalOperator(1/hitmap, broadcast='rightward')

    solution = pcg(A, b, M=M, disp=True, tol=1e-4, maxiter=200)
    #solution = pcg(A, b, disp=True, tol=1e-4, maxiter=200)
    x = pack.T*solution['x']
    x = solution['x']
    x[hitmap == 0] = np.nan
    sky_map = x.T

    return sky_map, hitmap


def get_signal(settings, bolo_name, segment):
    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.scanning_time_stamp)
    file_suffix = str(segment+1).zfill(4) + '.npy'
    ts = np.load(os.path.join(out_dir, bolo_name, 'ts_'+file_suffix))
    v = np.load(os.path.join(out_dir, bolo_name, 'vec_'+file_suffix))
    pol = np.load(os.path.join(out_dir, bolo_name, 'pol_'+file_suffix))
    return ts, v, pol

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
            signal, v, pol_ang= get_signal(settings, bolo_name, segment)
            if count%size is rank:
                print "Doing Bolo : ", bolo_name, "Segment : ", segment, "Rank : ", rank, "Count : ", count 
                sky_map, hitmap = make_map_from_signal(signal, v, pol_ang, bolo_name, segment)
            count+=1

    if rank is 0:
        write_map(settings, sky_map, hitmap)
        
if __name__=="__main__":
    from custom_settings import settings
    run_mpi(settings)
