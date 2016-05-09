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
    inv_cov_matrix[..., 1, 2] = sin_4 
    inv_cov_matrix[..., 2, 2] = 0.5*(n - cos_4) 
    inv_cov_matrix[..., 1, 0] = inv_cov_matrix[..., 0, 1]
    inv_cov_matrix[..., 2, 0] = inv_cov_matrix[..., 0, 2]
    inv_cov_matrix[..., 2, 1] = inv_cov_matrix[..., 1, 2]

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

    inv_cov_matrix_total *= 0.25

    hitmap = 4*inv_cov_matrix_total[..., 0, 0]

    bad_pix = hitmap<3

    inv_cov_matrix_total[bad_pix] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    cov_matrix = np.linalg.inv(inv_cov_matrix_total)

    sky_rec = np.sum(cov_matrix*b_total[..., None], axis=1).T

    sky_rec[..., bad_pix] = np.nan

    write_map(sky_rec, hitmap)


def get_local_bolo_segment_list(rank, num_processes):
    num_bolos = len(map_making_params.bolo_names)
    num_segments_per_bolo = map_making_params.segment_list.size
    num_total_segments = num_bolos*num_segments_per_bolo

    if num_total_segments % num_processes != 0:
        num_segments_per_process = num_total_segments/num_processes + 1
    else:
        num_segments_per_process = num_total_segments/num_processes
    
    bolo_list = []
    for bolo_name in map_making_params.bolo_names:
        bolo_list.extend([bolo_name]*num_segments_per_bolo)

    segment_list = np.arange(num_total_segments) % num_segments_per_bolo

    local_bolo_list = bolo_list[rank*num_segments_per_process : (rank + 1)*num_segments_per_process]
    local_segment_list = segment_list[rank*num_segments_per_process : (rank + 1)*num_segments_per_process]

    return local_bolo_list, local_segment_list, num_segments_per_process


def run_mpi():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size() 

    start_time_total = time.time()

    data_dir = os.path.join(map_making_params.global_output_dir, "scanning", map_making_params.scanning_time_stamp)

    npix = hp.nside2npix(map_making_params.nside_out)

    if rank == 0:
        inv_cov_matrix = np.zeros((npix, 3, 3))
        b_matrix = np.zeros((npix, 3))
    else:
        inv_cov_matrix = None
        b_matrix = None

    inv_cov_matrix_local = np.zeros((npix, 3, 3))
    b_matrix_local = np.zeros((npix, 3))

    local_bolo_list, local_segment_list, num_segments_per_process = get_local_bolo_segment_list(rank, size)
    print "Rank :", rank, "doing Bolos :", local_bolo_list, "and segments :", local_segment_list
    sys.stdout.flush()

    start_time = time.time()

    for bolo_name, segment in zip(local_bolo_list, local_segment_list):
        signal, v, pol_ang = get_signal(data_dir, bolo_name, segment)
        hitpix = hp.vec2pix(map_making_params.nside_out, v[...,0], v[...,1], v[...,2])
        b_matrix_local += get_b_matrix(hitpix, pol_ang, signal)
        inv_cov_matrix_local += get_inv_cov_matrix(hitpix, pol_ang)

    end_time = time.time()

    comm.Reduce(b_matrix_local, b_matrix, MPI.SUM, 0)
    comm.Reduce(inv_cov_matrix_local, inv_cov_matrix, MPI.SUM, 0)

    if rank == 0:

        del b_matrix_local
        del inv_cov_matrix_local

        inv_cov_matrix *= 0.25

        res_100 = [5748339, 5749950, 5749968, 5750147, 5751494, 5751521, 5751853, 5752077, 5754245, 5754291]
        res_01 = [16178, 18790, 19212, 19213, 19605, 19606, 19607, 20003, 20776, 20806]
        res_0 = [0, 1, 2, 3, 4 , 5, 6, 7, 8, 9, 10]

        hitmap = 4*inv_cov_matrix[..., 0, 0]

        bad_pix = hitmap<3

        inv_cov_matrix[bad_pix] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

        np.save("bad_pix_inv_100", inv_cov_matrix[res_100])
        np.save("bad_pix_inv_01", inv_cov_matrix[res_01])
        np.save("bad_pix_inv_0", inv_cov_matrix[res_0])
        np.save("bad_pix_b_100", b_matrix[res_100])
        np.save("bad_pix_b_01", b_matrix[res_01])
        np.save("bad_pix_b_0", b_matrix[res_0])

        cov_matrix = np.linalg.inv(inv_cov_matrix)

        np.save("bad_pix_100", cov_matrix[res_100])
        np.save("bad_pix_01", cov_matrix[res_01])
        np.save("bad_pix_0", cov_matrix[res_0])

        sky_rec = np.sum(cov_matrix*b_matrix[..., None], axis=1).T

        sky_rec[..., bad_pix] = np.nan

        write_map(sky_rec, hitmap)

        end_time_total = time.time()

        print "Average time taken per segment :", (end_time - start_time) / num_segments_per_process, "s"
        print "Total time taken :", (end_time_total - start_time_total), "s"


if __name__=="__main__":
    from custom_params import map_making_params
    action = sys.argv[1]

    if action=='run_mpi':
        from mpi4py import MPI
        run_mpi()

    if action=='run_serial':
        run_serial()
