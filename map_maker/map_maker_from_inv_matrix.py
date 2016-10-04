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
from simulation.lib.utilities.time_util import get_time_stamp 


def write_maps_and_config(sky_map, hitmap, bad_pix, recon_dir):
    hp.write_map(os.path.join(recon_dir, "reconstructed_map.fits"), sky_map)
    hitmap[bad_pix] = 0
    hp.write_map(os.path.join(recon_dir, "hitmap.fits"), hitmap)


def write_covariance_maps(maps, map_type, recon_dir):
    out_dir = os.path.join(recon_dir, map_type)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    if map_type == "partial_covariance_maps":
        map_legends = {"QQ" : (0,0), "QU" : (0,1), "UU" : (1,1)}
    else:
        map_legends = {"TT" : (0,0), "TQ" : (0,1), "TU" : (0,2), "QQ" : (1,1), "QU" : (1,2), "UU" : (2,2)}

    for leg in map_legends.keys():
        hp.write_map(os.path.join(out_dir, "map_" + leg + ".fits"), maps[..., map_legends[leg][0], map_legends[leg][1]])


def run_serial():
    start = time.time()

    npix = hp.nside2npix(config.nside_out)
    
    sim_dir = os.path.join(config.general_data_dir, config.sim_tag)
    recon_dir = os.path.join(config.general_data_dir, config.sim_tag, config.map_making_tag)

    inv_cov_matrix = np.empty((npix, 3, 3)) 
    b_matrix = np.load(os.path.join(recon_dir, "b_matrix.npy")) 
    inv_cov_matrix[..., 0, 0] = hp.read_map(os.path.join(recon_dir, "inverse_covariance_maps", "map_TT.fits"))
    inv_cov_matrix[..., 0, 1] = hp.read_map(os.path.join(recon_dir, "inverse_covariance_maps", "map_TQ.fits"))
    inv_cov_matrix[..., 0, 2] = hp.read_map(os.path.join(recon_dir, "inverse_covariance_maps", "map_TU.fits"))
    inv_cov_matrix[..., 1, 0] = inv_cov_matrix[..., 0, 1]
    inv_cov_matrix[..., 1, 1] = hp.read_map(os.path.join(recon_dir, "inverse_covariance_maps", "map_QQ.fits"))
    inv_cov_matrix[..., 1, 2] = hp.read_map(os.path.join(recon_dir, "inverse_covariance_maps", "map_QU.fits"))
    inv_cov_matrix[..., 2, 0] = inv_cov_matrix[..., 0, 2]
    inv_cov_matrix[..., 2, 1] = inv_cov_matrix[..., 1, 2]
    inv_cov_matrix[..., 2, 2] = hp.read_map(os.path.join(recon_dir, "inverse_covariance_maps", "map_UU.fits"))

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
    
    stop = time.time()

    prompter.prompt("Total time taken : %d" % (stop - start))


if __name__=="__main__":
    config_file = sys.argv[1]

    config = importlib.import_module("simulation.map_maker.config_files." + config_file).config

    run_serial()
