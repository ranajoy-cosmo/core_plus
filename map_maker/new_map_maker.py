#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import shutil
import time
from memory_profiler import profile
from mpi4py import MPI
from pysimulators.sparse import FSRBlockMatrix
from pysimulators import ProjectionOperator
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.timestream_simulation.new_sim_timestream import Bolo
from simulation.lib.utilities.prompter import prompt

#@profile
def get_signal(out_dir, bolo_name, segment):
    segment_name = str(segment+1).zfill(4)
    data_dir = os.path.join(out_dir, bolo_name, segment_name)
    ts = np.load(os.path.join(data_dir, "ts_signal.npy"))
    v = np.load(os.path.join(data_dir, "vector.npy"))
    pol = np.load(os.path.join(data_dir, "pol_ang.npy"))
    return ts, v, pol

#@profile
def get_inv_cov_matrix(hitpix, pol):
    nsamples = hitpix.size
    npix = hp.nside2npix(map_making_params.nside_out)

    n = np.bincount(hitpix, minlength=npix)
    cos_2 = np.bincount(hitpix, weights=np.cos(2*pol), minlength=npix)
    sin_2 = np.bincount(hitpix, weights=np.sin(2*pol), minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)
    sin_4 = np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)

    inv_cov_matrix = np.empty((npix, 3, 3))
    inv_cov_matrix[..., 0, 0] = n
    inv_cov_matrix[..., 0, 1] = cos_2 
    inv_cov_matrix[..., 0, 2] = sin_2 
    inv_cov_matrix[..., 1, 1] = 0.5*(n + cos_4) 
    inv_cov_matrix[..., 1, 2] = 0.5*sin_4 
    inv_cov_matrix[..., 2, 2] = 0.5*(n - cos_4) 
    inv_cov_matrix[..., 1, 0] = inv_cov_matrix[..., 0, 1]
    inv_cov_matrix[..., 2, 0] = inv_cov_matrix[..., 0, 2]
    inv_cov_matrix[..., 2, 1] = inv_cov_matrix[..., 1, 2]

    inv_cov_matrix *= 0.25

    return inv_cov_matrix

#@profile
def get_b_matrix(hitpix, pol, ts):
    nsamples = hitpix.size
    npix = hp.nside2npix(map_making_params.nside_out)
    
    matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.index[:, 0] = hitpix
    matrix.data.value[:, 0, 0, 0] = 0.5
    matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol)
    matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol)

    P = ProjectionOperator(matrix)

    b = P.T(ts)

    return b

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
def run_serial():

    data_dir = os.path.join(map_making_params.global_output_dir, "scanning", map_making_params.scanning_time_stamp)

    npix = hp.nside2npix(map_making_params.nside_out)
    inv_cov_matrix_total = np.zeros((npix, 3, 3))
    b_total = np.zeros((npix, 3))

    for bolo_name in map_making_params.bolo_names:
        for segment in map_making_params.segment_list:
            print "Doing bolo :", bolo_name, "and segment :", segment
            sys.stdout.flush()
            signal, v, pol_ang = get_signal(data_dir, bolo_name, segment)
            hitpix = hp.vec2pix(map_making_params.nside_out, v[...,0], v[...,1], v[...,2])
            start_time = time.time()
            b_total += get_b_matrix(hitpix, pol_ang, signal)
            end_time = time.time()
            #print "Time taken to get b matrix :", (end_time - start_time), "s"
            start_time = time.time()
            inv_cov_matrix_total += get_inv_cov_matrix(hitpix, pol_ang)
            end_time = time.time()
            #print "Time taken to get inverse covariance matrix :", (end_time - start_time), "s"
            sys.stdout.flush()

    hitmap = 4*inv_cov_matrix_total[..., 0, 0]

    bad_pix = hitmap<3

    inv_cov_matrix_total[bad_pix] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    cov_matrix = np.linalg.inv(inv_cov_matrix_total)

    sky_rec = np.sum(cov_matrix*b_total[..., None], axis=1).T

    sky_rec[..., bad_pix] = np.nan

    write_map(sky_rec, hitmap)



def run_mpi():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size() 

    #data_dir = os.path.join(map_making_params.global_output_dir, "scanning", map_making_params.scanning_time_stamp)

    npix = hp.nside2npix(config.nside_out)

    inv_cov_matrix_local = np.zeros((npix, 3, 3))
    b_matrix_local = np.zeros((npix, 3))

    bolo_segment_dict = get_local_bolo_segment_list(rank, size, map_making_params.bolo_list, map_making_params.segment_list)
    sys.stdout.flush()


    for bolo_name in bolo_segment_dict.keys():
        bolo = Bolo(bolo_name, config)
        for segment in bolo_segment_dict[bolo_name]:
            prompt("Rank : %d doing Bolo : %s and segment : %d" % (rank, bolo_name, segment))
            if config.simulate_ts:
                signal, v, pol_ang = bolo.simulate_timestream(...)
            else 
                signal, v, pol_ang = bolo.load_signal()
            hitpix = hp.vec2pix(map_making_params.nside_out, v[...,0], v[...,1], v[...,2])
            del v
            b_matrix_local += get_b_matrix(hitpix, pol_ang, signal)
            inv_cov_matrix_local += get_inv_cov_matrix(hitpix, pol_ang)
            time_out = time.time()


    if rank == 0:
        inv_cov_matrix = np.zeros((npix, 3, 3))
        b_matrix = np.zeros((npix, 3))
    else:
        inv_cov_matrix = None
        b_matrix = None

    comm.Reduce(b_matrix_local, b_matrix, MPI.SUM, 0)
    comm.Reduce(inv_cov_matrix_local, inv_cov_matrix, MPI.SUM, 0)

    if rank == 0:

        del b_matrix_local
        del inv_cov_matrix_local

        inv_cov_matrix *= 0.25

        hitmap = 4*inv_cov_matrix[..., 0, 0]

        bad_pix = hitmap<3

        inv_cov_matrix[bad_pix] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

        cov_matrix = np.linalg.inv(inv_cov_matrix)

        sky_rec = np.sum(cov_matrix*b_matrix[..., None], axis=1).T

        sky_rec[..., bad_pix] = np.nan

        write_map(sky_rec, hitmap)


if __name__=="__main__":
    config_file = sys.argv[1]

    config = importlib.import_module(config_file).config

    if config.run_type=='run_mpi':
        from mpi4py import MPI
        run_mpi()

    if config.run_type=='run_serial':
        run_serial()
