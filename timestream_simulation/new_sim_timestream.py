#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
import copy
import os
import shutil
import importlib
import time
from memory_profiler import profile
from simulation.lib.quaternion import quaternion
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
from simulation.beam.beam_kernel import get_beam, display_beam_settings
from simulation.lib.utilities.time_util import get_time_stamp
import simulation.lib.numericals.filters as filters
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.lib.utilities.generic_class import Generic
from simulation.beam.new_beam_kernel import Beam
import simulation.lib.utilities.prompter as prompter


def display_params():
    display_string = ""
    display_string += "Alpha : %f degrees\n" % (config.alpha)
    display_string += "Beta : %f degrees\n" % (config.beta)
    t_flight = config.t_segment*len(config.segment_list)
    display_string += "T flight : %f hours / %f days" % (t_flight/60.0/60.0, t_flight/60.0/60.0/24.0)
    display_string += "T segment : %f hours / %f days" % (config.t_segment/60.0/60.0, config.t_segment/60.0/60.0/24)
    display_string += "T precession : %f hours" % (config.t_prec/60.0/60.0)
    display_string += "T spin : %f seconds" % (config.t_spin)
    display_string += "Scan sampling rate : %f Hz" % (config.sampling_rate)
    display_string += "Theta co : %f arcmin" % (config.theta_co)
    display_string += "Theta cross : %f arcmin" % (config.theta_cross)
    display_string += "Oversampling rate : %d" % (config.oversampling_rate)
    display_string += "Scan resolution for beam integration : %f arcmin" % (config.scan_resolution)
    display_string += "Pixel size for NSIDE = %d : %f arcmin" % (config.nside, hp.nside2resol(config.nside, arcmin=True))
    n_steps = int(config.t_segment*config.sampling_rate)*config.oversampling_rate 
    display_string += "#Samples per segment : %d", %(n_steps)

    prompt(display_string, False if not config.action=="display_params" else True)


def make_data_dirs(data_dir):
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    info_dir = os.path.join(data_dir, "run_info")
    if not os.path.exists(info_dir):
        os.makedirs(info_dir)
    shutil.copy("default_config.py", info_dir)
    shutil.copy(config_file + ".py", info_dir)
    shutil.copy("new_sim_timestream.py", info_dir)
    shutil.copy(config.bolo_config_file + ".py", info_dir)
    for bolo in config.bolo_list:
        bolo_dir = os.path.join(data_dir, bolo)
        if not os.path.exists(bolo_dir):
            os.makedirs(bolo_dir)
        for segment in config.segment_list:
            segment_name = str(segment+1).zfill(4)
            segment_dir = os.path.join(data_dir, bolo, segment_name)
            if not os.path.exists(segment_dir):
                os.makedirs(segment_dir)
            
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* The master section which distributes the data and collects it
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*


def run_mpi():
    sim_tag = comm.bcast(config.sim_tag, root=0)         #All processes get the tag of Rank 0. Important when the tag is a time stamp
    scan_tag = comm.bcast(config.scan_tag, root=0)         #All processes get the tag of Rank 0. Important when the tag is a time stamp

    data_dir = check_and_make_data_dirs(tag)

    comm.Barrier()

    npix_out = hp.nside2npix(config.nside_out)
    inv_cov_matrix_local = np.zeros((6, npix_out))
    b_matrix_local = np.zeros((npix_out, 3))

    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)

    for bolo_name in bolo_segment_dict.keys():
        bolo = Bolo(bolo_name, config, data_dir)
        for segment in bolo_segment_dict[bolo_name]:
            prompt("Doing Bolo : %s Segment : %d Rank : %d" % (bolo_name, segment+1, rank))
            if config.do_pencil_beam:
                signal, v, pol_ang = bolo.simulate_timestream(segment)
            else:
                signal, v, pol_ang = bolo.simulate_timestream_beamed(segment)
            hitpix = hp.vec2pix(config.nside_out, v[...,0], v[...,1], v[...,2])
            del v
            b_matrix_local += get_b_matrix(hitpix, pol_ang, signal)
            inv_cov_matrix_local += get_inv_cov_matrix(hitpix, pol_ang)

    del hitpix, pol_ang, signal
    
    if config.make_scanned_map:
        if rank == 0:
            hitmap = np.zeros(hp.nside2npix(config.nside), dtype=np.float32)
        else:
            hitmap = None

        comm.Reduce(hitmap_local, hitmap, MPI.SUM, 0)
        del hitmap_local

        if rank == 0:
            valid = hitmap>0
            if config.do_only_T:
                bolo.sky_map[~valid] = np.nan
            else:
                bolo.sky_map[...,~valid] = np.nan
            hp.write_map(os.path.join(data_dir, "scanned_map.fits"), bolo.sky_map)
            hp.write_map(os.path.join(data_dir, "hitmap_in.fits"), hitmap)


def run_check(config, verbose=True):
    if config.noise_only:
        if not config.do_pencil_beam:
            if verbose:
                prompter.prompt_warning("No beam convolution required for noise only simulation.")
            config.do_pencil_beam = True

    if config.do_pencil_beam:
        if config.oversampling_rate != 1:
            if verbose:
                prompter.prompt_warning("No oversampling required for pencil beam")
            config.oversampling_rate = 1


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Main function definition. This is where the code begins when executed
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    config_file = sys.argv[1]
    run_type = sys.argv[2]

    config = importlib.import_module(config_file).config

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
