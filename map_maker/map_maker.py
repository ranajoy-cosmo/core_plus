#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import shutil
import time
import importlib
#from memory_profiler import profile
from pysimulators.sparse import FSRBlockMatrix
from pysimulators import ProjectionOperator
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.timestream_simulation.bolo import Bolo
import simulation.lib.utilities.prompter as prompter
from simulation.lib.utilities.time_util import get_time_stamp 

def get_inv_cov_matrix(v, pol, ts, inv_cov_matrix, b_matrix):
    nsamples = hitpix.size
    npix = hp.nside2npix(config.nside_out)

    hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])

    n = np.bincount(hitpix, minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)
    cos_2_pol = np.cos(2*pol)
    sin_2_pol = np.sin(2*pol)

    if config.pol_type == "TQU"
        inv_cov_matrix[..., 0] += 0.25*n
        inv_cov_matrix[..., 1] += 0.25*np.bincount(hitpix, weights=cos_2_pol, minlength=npix)
        inv_cov_matrix[..., 2] += 0.25*np.bincount(hitpix, weights=sin_2_pol, minlength=npix)
        inv_cov_matrix[..., 3] += 0.25*0.5*(n + cos_4) 
        inv_cov_matrix[..., 4] += 0.25*0.5*np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)
        inv_cov_matrix[..., 5] += 0.25*0.5*(n - cos_4) 

        b_matrix[..., 0] += np.bincount(hitpix, 0.5*ts, minlength=npix)
        b_matrix[..., 1] += np.bincount(hitpix, 0.5*cos_2_pol*ts, minlength=npix)
        b_matrix[..., 2] += np.bincount(hitpix, 0.5*sin_2_pol*ts, minlength=npix)
    else:
        inv_cov_matrix[..., 0] += 0.25*n
        inv_cov_matrix[..., 1] += 0.25*0.5*(n + cos_4) 
        inv_cov_matrix[..., 2] += 0.25*0.5*np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)
        inv_cov_matrix[..., 3] += 0.25*0.5*(n - cos_4) 

        b_matrix[..., 0] += np.bincount(hitpix, 0.5*cos_2_pol*ts, minlength=npix)
        b_matrix[..., 1] += np.bincount(hitpix, 0.5*sin_2_pol*ts, minlength=npix)


def write_maps_and_config(sky_map, hitmap, bad_pix, recon_dir):
    hp.write_map(os.path.join(recon_dir, "reconstructed_map.fits"), sky_map)
    hitmap[bad_pix] = 0
    hp.write_map(os.path.join(recon_dir, "hitmap.fits"), hitmap)


def write_covariance_maps(maps, map_type, recon_dir):
    out_dir = os.path.join(recon_dir, map_type)
    os.makedirs(out_dir)

    if map_type == "partial_covariance_maps" || config.pol_type == "pol_only":
        map_legends = {"QQ" : 0, "QU" : 1, "UU" : 2}
    else:
        map_legends = {"TT" : 0, "TQ" : 1, "TU" : 2, "QQ" : 3, "QU" : 4, "UU" : 5}

    for leg in map_legends.keys():
        hp.write_map(os.path.join(out_dir, "map_" + leg + ".fits"), maps[..., map_legends[leg][0], map_legends[leg][1]])
        #np.save(os.path.join(out_dir, "map_" + leg), maps[..., map_legends[leg][0], map_legends[leg][1]])


def get_signal(segment, bolo_a, bolo_b):
    if config.mode == "single_bolo"
        if config.simulate_ts:
            signal, v, pol_ang = bolo_a.simulate_timestream(segment)
        else:
            signal, v, pol_ang = bolo_a.read_timestream(segment)
    else:
        if config.simulate_ts:
            signal, v, pol_ang = bolo_a.simulate_timestream(segment)
            signal -= bolo_b.simulate_timestream(segment, return_field=["signal"])
        else:
            signal, v, pol_ang = bolo_a.read_timestream(segment)
            signal -= bolo_b.read_timestream(segment, return_field=["signal"])

    return signal, v, pol_ang


def run_mpi():
    start = time.time()
    npix = hp.nside2npix(config.nside_out)

    if config.pol_type == "TQU"
        inv_cov_matrix_local = np.zeros((npix, 6))
        b_matrix_local = np.zeros((npix, 3))
    else:
        inv_cov_matrix_local = np.zeros((npix, 3))
        b_matrix_local = np.zeros((npix, 2))

    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)

    print "Rank :", rank, ", Bolos and Segments :", bolo_segment_dict
    comm.Barrier()

    if rank == 0:
        recon_dir = make_data_dirs() 

    for bolo_name in bolo_segment_dict.keys():
        if config.mode = "single_bolo":
            bolo_1 = Bolo(bolo_name, config)
            bolo_2 = None
        else:
            bolo_1 = Bolo(bolo_name + 'a', config)
            bolo_2 = Bolo(bolo_name + 'b', config)
        for segment in bolo_segment_dict[bolo_name]:
            segment_start = time.time()
            prompter.prompt("Rank : %d doing Bolo : %s and segment : %d" % (rank, bolo_name, segment))
            signal, v, pol_ang = get_signal(segment, bolo_1, bolo_2)
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            get_inv_cov_matrix(hitpix, pol_ang, signal, inv_cov_matrix_local, b_matrix_local)
            segment_stop = time.time()
            prompter.prompt("Rank : %d doing Bolo : %s and segment : %d and time taken : %d" % (rank, bolo_name, segment, segment_stop - segment_start))

    if config.pol_type == "TQU"
        inv_cov_matrix= np.zeros((npix, 6))
        b_matrix= np.zeros((npix, 3))
    else:
        inv_cov_matrix= np.zeros((npix, 3))
        b_matrix= np.zeros((npix, 2))

    comm.Reduce(b_matrix_local, b_matrix, MPI.SUM, 0)
    comm.Reduce(inv_cov_matrix_local, inv_cov_matrix, MPI.SUM, 0)

    del b_matrix_local
    del inv_cov_matrix_local

    hitpix = 4*inv_cov_matrix[..., 0]

    if config.pol_type == "TQU"
        inv_cov_matrix[hitpix<3] = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 1.0])
    else:
        inv_cov_matrix[hitpix<3] = np.array([1.0, 0.0, 1.0])
    write_covariance_maps(inv_cov_matrix, "inverse_covariance_maps", recon_dir)


        cov_matrix = np.linalg.inv(inv_cov_matrix)
        write_covariance_maps(cov_matrix, "covariance_maps", recon_dir)

        if "partial_covariance_maps" in config.map_making_data_products:
            cov_matrix_partial = np.linalg.inv(inv_cov_matrix[..., 1:, 1:])
            write_covariance_maps(cov_matrix_partial, "partial_covariance_maps", recon_dir)

        sky_rec = np.sum(cov_matrix*b_matrix[..., None], axis=1).T

        sky_rec[..., bad_pix] = np.nan

        write_maps_and_config(sky_rec, 4*inv_cov_matrix[..., 0, 0], bad_pix, recon_dir)
        
        stop = time.time()

        prompter.prompt("Total time taken : %d" % (stop - start))


def run_serial():
    start = time.time()
    npix = hp.nside2npix(config.nside_out)

    if config.pol_type == "TQU"
        inv_cov_matrix_unravel = np.zeros((npix, 6))
        b_matrix = np.zeros((npix, 3))
    else:
        inv_cov_matrix_unravel = np.zeros((npix, 3))
        b_matrix = np.zeros((npix, 2))

    recon_dir = make_data_dirs() 

    for bolo_name in config.bolo_list: 
        if config.mode = "single_bolo":
            bolo_1 = Bolo(bolo_name, config)
            bolo_2 = None
        else:
            bolo_1 = Bolo(bolo_name + 'a', config)
            bolo_2 = Bolo(bolo_name + 'b', config)
        for segment in config.segment_list: 
            segment_start = time.time()
            prompter.prompt("Bolo : %s and segment : %d" % (bolo_name, segment))
            signal, v, pol_ang = get_signal(segment, bolo_1, bolo_2)
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            get_inv_cov_matrix(hitpix, pol_ang, signal, inv_cov_matrix_unravel, b_matrix)
            segment_stop = time.time()
            prompter.prompt("Bolo : %s and segment : %d and time taken : %d" % (bolo_name, segment, segment_stop - segment_start))

    if config.pol_type == "TQU"
        inv_cov_matrix= np.zeros((npix, 3, 3))
    else:
        inv_cov_matrix= np.zeros((npix, 3))

    hitpix = 4*inv_cov_matrix[..., 0]

    if config.pol_type == "TQU"
        inv_cov_matrix[hitpix<3] = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 1.0])
    else:
        inv_cov_matrix[hitpix<3] = np.array([1.0, 0.0, 1.0])
    write_covariance_maps(inv_cov_matrix, "inverse_covariance_maps", recon_dir)


        cov_matrix = np.linalg.inv(inv_cov_matrix)
        write_covariance_maps(cov_matrix, "covariance_maps", recon_dir)

        if "partial_covariance_maps" in config.map_making_data_products:
            cov_matrix_partial = np.linalg.inv(inv_cov_matrix[..., 1:, 1:])
            write_covariance_maps(cov_matrix_partial, "partial_covariance_maps", recon_dir)

        sky_rec = np.sum(cov_matrix*b_matrix[..., None], axis=1).T

        sky_rec[..., bad_pix] = np.nan

        write_maps_and_config(sky_rec, 4*inv_cov_matrix[..., 0, 0], bad_pix, recon_dir)
        
        stop = time.time()

        prompter.prompt("Total time taken : %d" % (stop - start))


def make_data_dirs():
    sim_dir = os.path.join(config.general_data_dir, config.sim_tag)
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    time_stamp = get_time_stamp()

    recon_dir = os.path.join(config.general_data_dir, config.sim_tag, config.map_making_tag)
    config_dir = os.path.join(recon_dir, "config_files")
    default_config_file = os.path.join(config.base_dir, "map_maker/config_files/default_config.py")
    current_config_file = os.path.join(config.base_dir, "map_maker/config_files", config_file + ".py")

    if os.path.exists(recon_dir):
        if config.map_making_action == "new":
            shutil.rmtree(recon_dir)
            os.makedirs(recon_dir)
            os.makedirs(config_dir)
        else:
            pass
    else:
        os.makedirs(recon_dir)
        os.makedirs(config_dir)

    shutil.copy(default_config_file, config_dir)
    shutil.copyfile(current_config_file, os.path.join(config_dir, config_file + time_stamp + ".py"))
    
    if config.timestream_data_products and config.simulate_ts:
        scan_dir = os.path.join(sim_dir, config.scan_tag)
        if not os.path.exists(scan_dir):
            os.makedirs(scan_dir)
        config_dir = os.path.join(scan_dir, "config_files")
        if not os.path.exists(config_dir):
            os.makedirs(config_dir)
        default_config_file = "/global/homes/b/banerji/simulation/map_maker/config_files/default_config.py"
        shutil.copy(default_config_file, config_dir)
        current_config_file = os.path.join("/global/homes/b/banerji/simulation/map_maker/config_files", config_file + '.py') 
        shutil.copy(current_config_file, config_dir)


    return recon_dir


if __name__=="__main__":
    config_file = sys.argv[1]
    run_type = sys.argv[2]

    config = importlib.import_module("simulation.map_maker.config_files." + config_file).config

    if run_type=='run_mpi':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        run_mpi()

    if run_type=='run_serial':
        run_serial()
