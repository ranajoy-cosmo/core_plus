#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import shutil
import time
import importlib
#from memory_profiler import profile
from mpi4py import MPI
from pysimulators.sparse import FSRBlockMatrix
from pysimulators import ProjectionOperator
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.timestream_simulation.bolo import Bolo
import simulation.lib.utilities.prompter as prompter
from simulation.lib.utilities.time_util import get_time_stamp 

def get_inv_cov_matrix(hitpix, pol):
    nsamples = hitpix.size
    npix = hp.nside2npix(config.nside_out)

    n = np.bincount(hitpix, minlength=npix)
    cos_2 = np.bincount(hitpix, weights=np.cos(2*pol), minlength=npix)
    sin_2 = np.bincount(hitpix, weights=np.sin(2*pol), minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)
    sin_4 = np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)

    inv_cov_matrix = np.empty((npix, 2, 2))
    inv_cov_matrix[..., 0, 0] = 0.5*(n + cos_4) 
    inv_cov_matrix[..., 0, 1] = 0.5*sin_4 
    inv_cov_matrix[..., 1, 0] = inv_cov_matrix[..., 0, 1]
    inv_cov_matrix[..., 1, 1] = 0.5*(n - cos_4) 

    return 0.25*inv_cov_matrix


def get_b_matrix(hitpix, pol, ts):
    nsamples = hitpix.size
    npix = hp.nside2npix(config.nside_out)
    
    matrix = FSRBlockMatrix((nsamples, npix*2), (1, 2), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.index[:, 0] = hitpix
    matrix.data.value[:, 0, 0, 0] = 0.5*np.cos(2*pol)
    matrix.data.value[:, 0, 0, 1] = 0.5*np.sin(2*pol)

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

    map_legends = {"QQ" : (0,0), "QU" : (0,1), "UU" : (1,1)}

    for leg in map_legends.keys():
        hp.write_map(os.path.join(out_dir, "map_" + leg), maps[..., map_legends[leg][0], map_legends[leg][1]])
        #np.save(os.path.join(out_dir, "map_" + leg), maps[..., map_legends[leg][0], map_legends[leg][1]])


def run_mpi():
    start = time.time()
    npix = hp.nside2npix(config.nside_out)

    inv_cov_matrix_local = np.zeros((npix, 2, 2))
    b_matrix_local = np.zeros((npix, 2))
    hitmap_local = np.zeros(npix)

    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)

    if rank == 0:
        recon_dir = make_data_dirs() 

    comm.Barrier()

    #alpha = dict(np.load("/global/homes/h/hoangapc/simulation/output/100_bolo_run/20_bolos_scan_with_noise/alpha.npy"))
    #alpha = {'bolo_0010': 3.3436891282375049e-06, 'bolo_0008': 6.9888745498721233e-05, 'bolo_0009': 9.1838785556966411e-05, 'bolo_0001': 8.4746935102100737e-05, 'bolo_0002': 2.5493732551259379e-06, 'bolo_0003': 3.3242380731229615e-05, 'bolo_0004': 1.0465765890885307e-05, 'bolo_0005': 0.00019143450831889689, 'bolo_0006': -1.2783283358750637e-05, 'bolo_0007': 6.1247427220433365e-05}
    #alpha = {'bolo_0010': 1.6718445641187505e-06, 'bolo_0008': 3.4944372749360596e-05, 'bolo_0009': 4.5919392778483212e-05, 'bolo_0001': 4.2373467551050348e-05, 'bolo_0002': 1.2746866275629672e-06, 'bolo_0003': 1.6621190365614804e-05, 'bolo_0004': 5.2328829454426507e-06, 'bolo_0005': 9.571725415944847e-05, 'bolo_0006': -6.3916416793753186e-06, 'bolo_0007': 3.0623713610216682e-05}
    alpha = {'bolo_0010': 1.7572309733761022e-06, 'bolo_0008': 3.4952490679335271e-05, 'bolo_0009': 4.5895232229898289e-05, 'bolo_0001': 4.2347476167810334e-05, 'bolo_0002': 1.2964126315678174e-06, 'bolo_0003': 1.6557738825949474e-05, 'bolo_0004': 5.2219631868142547e-06, 'bolo_0005': 9.553259540007582e-05, 'bolo_0006': -6.3018957304171546e-06, 'bolo_0007': 3.0564560757508306e-05}
    bolo_350 = Bolo("bolo_00350", config)
    for bolo_name in bolo_segment_dict.keys():
        bolo_a_name = bolo_name + 'a'
        bolo_b_name = bolo_name + 'b'
        bolo_a = Bolo(bolo_a_name, config)
        bolo_b = Bolo(bolo_b_name, config)
        for segment in bolo_segment_dict[bolo_name]:
            prompter.prompt("Rank : %d doing Bolo : %s and segment : %d" % (rank, bolo_name, segment))
            signal_b, v, pol_ang = bolo_b.read_timestream(segment)
            signal_a, v, pol_ang = bolo_a.read_timestream(segment)
            signal_350, v, pol_ang = bolo_350.read_timestream(segment)
            signal_diff = 0.5*(signal_a - signal_b) - alpha[bolo_name]*signal_350 
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            b_matrix_local += get_b_matrix(hitpix, pol_ang, signal_diff)
            inv_cov_matrix_local += get_inv_cov_matrix(hitpix, pol_ang)
            hitmap_local += np.bincount(hitpix, minlength=npix)

    inv_cov_matrix = np.zeros((npix, 2, 2))
    b_matrix = np.zeros((npix, 2))
    hitmap = np.zeros(npix)

    comm.Reduce(b_matrix_local, b_matrix, MPI.SUM, 0)
    comm.Reduce(inv_cov_matrix_local, inv_cov_matrix, MPI.SUM, 0)
    comm.Reduce(hitmap_local, hitmap, MPI.SUM, 0)

    if rank == 0:
        del b_matrix_local
        del inv_cov_matrix_local

        bad_pix = hitmap<3

        inv_cov_matrix[bad_pix] = np.array([[1.0, 0.0], [0.0, 1.0]])
        write_covariance_maps(inv_cov_matrix, "inverse_covariance_maps", recon_dir)

        cov_matrix = np.linalg.inv(inv_cov_matrix)
        write_covariance_maps(cov_matrix, "partial_covariance_maps", recon_dir)

        """
        if "partial_covariance_maps" in config.map_making_data_products:
            cov_matrix_partial = np.linalg.inv(inv_cov_matrix[..., 1:, 1:])
            write_covariance_maps(cov_matrix_partial, "partial_covariance_maps", recon_dir)
        """

        sky_rec = np.sum(cov_matrix*b_matrix[..., None], axis=1).T

        #sky_rec[..., bad_pix] = np.nan

        write_maps_and_config(sky_rec, hitmap, bad_pix, recon_dir)
        
        stop = time.time()

        prompter.prompt("Total time taken : %d" % (stop - start))


def make_data_dirs():
    sim_dir = os.path.join(config.general_data_dir, config.sim_tag)
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    time_stamp = get_time_stamp()

    recon_dir = os.path.join(config.general_data_dir, config.sim_tag, config.map_making_tag)
    config_dir = os.path.join(recon_dir, "config_files")
    default_config_file = "/global/homes/b/banerji/simulation/map_maker/config_files/default_config.py"
    current_config_file = os.path.join("/global/homes/b/banerji/simulation/map_maker/config_files", config_file + ".py")

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

    #shutil.copy(default_config_file, config_dir)
    #shutil.copyfile(current_config_file, os.path.join(config_dir, config_file + time_stamp + ".py"))
    
    if config.timestream_data_products:
        scan_dir = os.path.join(sim_dir, config.scan_tag)
        if not os.path.exists(scan_dir):
            os.makedirs(scan_dir)
        config_dir = os.path.join(scan_dir, "config_files")
        if not os.path.exists(config_dir):
            os.makedirs(config_dir)
        default_config_file = "/global/homes/b/banerji/simulation/timestream_simulation/config_files/default_config.py"
        #shutil.copy(default_config_file, config_dir)
        current_config_file = os.path.join("/global/homes/b/banerji/simulation/timestream_simulation/config_files", config_file + '.py') 
        #shutil.copy(current_config_file, config_dir)


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
