#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
from mpi4py import MPI
from pysimulators.sparse import FSRBlockMatrix
from pysimulators import ProjectionOperator

def get_signal(out_dir, bolo_name, segment):
    segment_name = str(segment+1).zfill(4)
    data_dir = os.path.join(out_dir, bolo_name, segment_name)
    ts = np.load(os.path.join(data_dir, "ts_signal.npy"))
    v = np.load(os.path.join(data_dir, "vector.npy"))
    pol = np.load(os.path.join(data_dir, "pol_ang.npy"))
    return ts, v, pol

def get_inv_cov_matrix(hitpix, pol):
    nsamples = hit_pix.size
    npix = hp.nside2npix(map_making_params.nside_out)

    n = np.bincount(hit_pix, minlength=npix)
    cos_2 = np.bincount(hit_pix, weights=np.cos(2*pol), minlength=npix)
    sin_2 = np.bincount(hit_pix, weights=np.sin(2*pol), minlength=npix)
    cos_4 = np.bincount(hit_pix, weights=np.cos(4*pol), minlength=npix)
    sin_4 = np.bincount(hit_pix, weights=np.sin(4*pol), minlength=npix)

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

def get_b_matrix(hitpix, pol, ts):
    nsamples = hit_pix.size
    npix = hp.nside2npix(map_making_params.nside_out)
    
    matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.index[:, 0] = hitpix
    matrix.data.value[:, 0, 0, 0] = 0.5
    matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol)
    matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol)

    P = ProjectionOperator(matrix)

    b = P.T(ts)

    return b

def run_serial():

    data_dir = os.path.join(map_making_params.global_output_dir, "scanning", map_making_params.scanning_time_stamp)

    npix = hp.nside2npix(map_making_params.nside_out)
    inv_cov_matrix_total = np.zeros((npix, 3, 3))
    b_total = np.zeros((npix, 3))

    for bolo_name in map_making_params.bolo_names:
        for segment in map_making_params.segment_list:
            signal, v, pol_ang = get_signal(data_dir, bolo_name, segment)
            hitpix = hp.vec2pix(map_making_params.nside_out, v[...,0], v[...,1], v[...,2])
            del v
            b_total += get_b_matrix(hitpix, pol_ang, signal)
            inv_cov_matrix_total += get_inv_cov_matrix(hitpix, pol_ang)

    inv_cov_matrix *= 0.25

    del signal
    del pol_ang

    hitmap = 4*inv_cov_matrix_total[..., 0, 0]

    bad_pix = hitmap<3

    inv_cov_matrix_total[bad_pix] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    cov_matrix = np.linalg.inv(inv_cov_matrix_total)

    del inv_cov_matrix_total

    sky_rec = np.sum(cov_matrix*b_total[..., None], axis=1).T

    return hitmap, sky_rec

if __name__=="__main__":
    from custom_params import map_making_params
    hitmap, sky_rec = run_serial()
