#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import shutil
import time
import importlib
from memory_profiler import profile
from mpi4py import MPI
from pysimulators.sparse import FSRBlockMatrix
from pysimulators import ProjectionOperator
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.timestream_simulation.new_sim_timestream import Bolo
import simulation.lib.utilities.prompter as prompter

def get_inv_cov_matrix(hitpix, pol):
    nsamples = hitpix.size
    npix = hp.nside2npix(config.nside_out)

    n = np.bincount(hitpix, minlength=npix)
    cos_2 = np.bincount(hitpix, weights=np.cos(2*pol), minlength=npix)
    sin_2 = np.bincount(hitpix, weights=np.sin(2*pol), minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)
    sin_4 = np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)

    inv_cov_matrix = np.empty((npix, 3, 3))
    inv_cov_matrix[..., 0, 0] = n
    inv_cov_matrix[..., 0, 1] = cos_2 
    inv_cov_matrix[..., 0, 2] = sin_2 
    inv_cov_matrix[..., 1, 0] = inv_cov_matrix[..., 0, 1]
    inv_cov_matrix[..., 1, 1] = 0.5*(n + cos_4) 
    inv_cov_matrix[..., 1, 2] = 0.5*sin_4 
    inv_cov_matrix[..., 2, 0] = inv_cov_matrix[..., 0, 2]
    inv_cov_matrix[..., 2, 1] = inv_cov_matrix[..., 1, 2]
    inv_cov_matrix[..., 2, 2] = 0.5*(n - cos_4) 

    inv_cov_matrix *= 0.25

    return inv_cov_matrix

def get_b_matrix(hitpix, pol, ts):
    nsamples = hitpix.size
    npix = hp.nside2npix(config.nside_out)
    
    matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.index[:, 0] = hitpix
    matrix.data.value[:, 0, 0, 0] = 0.5
    matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol)
    matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol)

    P = ProjectionOperator(matrix)

    b = P.T(ts)

    return b

def write_maps_and_config(sky_map, hitmap, bad_pix, recon_dir):
    hp.write_map(os.path.join(recon_dir, "reconstructed_map.fits"), sky_map)
    hitmap[bad_pix] = 0
    hp.write_map(os.path.join(recon_dir, "hitmap.fits"), hitmap)


def write_covariance_maps(maps, map_type, recon_dir):
    out_dir = os.path.join(recon_dir, map_type)
    os.makedirs(out_dir)

    if map_type == "partial_covariance_maps":
        map_legends = {"EE" : (0,0), "EQ" : (0,1), "QQ" : (1,1)}
    else:
        map_legends = {"TT" : (0,0), "TE" : (0,1), "TQ" : (0,2), "EE" : (1,1), "EQ" : (1,2), "QQ" : (2,2)}

    for leg in map_legends.keys():
        np.save(os.path.join(out_dir, "map_" + leg), maps[..., map_legends[leg][0], map_legends[leg][1]])

def run_mpi():
    npix = hp.nside2npix(config.nside_out)

    inv_cov_matrix_local = np.zeros((npix, 3, 3))
    b_matrix_local = np.zeros((npix, 3))

    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)

    if rank == 0:
        recon_dir = make_data_dirs() 

    for bolo_name in bolo_segment_dict.keys():
        bolo = Bolo(bolo_name, config)
        for segment in bolo_segment_dict[bolo_name]:
            prompter.prompt("Rank : %d doing Bolo : %s and segment : %d" % (rank, bolo_name, segment))
            if config.simulate_ts:
                signal, v, pol_ang = bolo.simulate_timestream(segment)
            else:
                signal, v, pol_ang = bolo.read_timestream(segment)
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            b_matrix_local += get_b_matrix(hitpix, pol_ang, signal)
            inv_cov_matrix_local += get_inv_cov_matrix(hitpix, pol_ang)


    inv_cov_matrix = np.zeros((npix, 3, 3))
    b_matrix = np.zeros((npix, 3))

    comm.Reduce(b_matrix_local, b_matrix, MPI.SUM, 0)
    comm.Reduce(inv_cov_matrix_local, inv_cov_matrix, MPI.SUM, 0)

    if rank == 0:
        del b_matrix_local
        del inv_cov_matrix_local

        bad_pix = 4*inv_cov_matrix[..., 0, 0]<3

        inv_cov_matrix[bad_pix] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        write_covariance_maps(inv_cov_matrix, "inverse_covariance_maps", recon_dir)

        cov_matrix = np.linalg.inv(inv_cov_matrix)
        write_covariance_maps(cov_matrix, "covariance_maps", recon_dir)

        if "partial_covariance_maps" in config.map_making_data_products:
            cov_matrix_partial = np.linalg.inv(inv_cov_matrix[..., 1:, 1:])
            write_covariance_maps(cov_matrix_partial, "partial_covariance_maps", recon_dir)

        sky_rec = np.sum(cov_matrix*b_matrix[..., None], axis=1).T

        sky_rec[..., bad_pix] = np.nan

        write_maps_and_config(sky_rec, 4*inv_cov_matrix[..., 0, 0], bad_pix, recon_dir)


def make_data_dirs():
    sim_dir = os.path.join(config.general_data_dir, config.sim_tag)
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    recon_dir = os.path.join(config.general_data_dir, config.sim_tag, config.map_making_tag)
    if os.path.exists(recon_dir):
        shutil.rmtree(recon_dir)
    os.makedirs(recon_dir)
    config_dir = os.path.join(recon_dir, "config_files")
    os.makedirs(config_dir)
    default_config_file = "/global/homes/b/banerji/simulation/map_maker/config_files/default_config.py"
    shutil.copy(default_config_file, config_dir)
    current_config_file = os.path.join("/global/homes/b/banerji/simulation/map_maker/config_files", config_file + ".py")
    shutil.copy(current_config_file, config_dir)
    
    if config.timestream_data_products:
        scan_dir = os.path.join(sim_dir, config.scan_tag)
        if not os.path.exists(scan_dir):
            os.makedirs(scan_dir)
        config_dir = os.path.join(scan_dir, "config_files")
        if not os.path.exists(config_dir):
            os.makedirs(config_dir)
        default_config_file = "/global/homes/b/banerji/simulation/timestream_simulation/config_files/default_config.py"
        shutil.copy(default_config_file, config_dir)
        current_config_file = os.path.join("/global/homes/b/banerji/simulation/timestream_simulation/config_files", config_file + '.py') 
        shutil.copy(current_config_file, config_dir)


    return recon_dir

def get_covariance_matrix(inv_cov_matrix):
    nv_cov_matrix_local = None 


def get_local_pix_range():
    npix = hp.nside2npix(config.nside_out)
    u_npix_pp = npix/size
    add_npix_0 = npix%u_npix_pp

    start = rank*u_npix_pp
    stop = (rank + 1)*u_npix_pp + add_npix_0
    if rank != 0:
        start += add_npix_0
    
    return start, stop

if __name__=="__main__":
    config_file = sys.argv[1]
    run_type = sys.argv[2]

    config = importlib.import_module("simulation.map_maker.config_files." + config_file).config

    if run_type=='run_check':
        run_check()

    if run_type=='run_mpi':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        run_mpi()

    if run_type=='run_serial':
        run_serial()
