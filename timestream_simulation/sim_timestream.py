#!/usr/bin/env python

import numpy as np
import healpy as hp
import sys
import os
import shutil
import importlib
import time
from memory_profiler import profile
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
import simulation.lib.utilities.prompter as prompter
from simulation.timestream_simulation.bolo import Bolo


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* The master section which distributes the data and collects it
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*


def run_mpi():
    if rank == 0:
        make_data_dirs()

    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)
    comm.Barrier()
    time.sleep(0.1*rank)
    print "Rank :", rank, "Local bolo segment list :\n", bolo_segment_dict
    comm.Barrier()

    tot_seg = 0
    for keys in bolo_segment_dict.keys():
        tot_seg += len(bolo_segment_dict[keys])

    if "hitmap" in config.timestream_data_products:
        hitmap_local = np.zeros(hp.nside2npix(config.nside_in))

    for bolo_name in bolo_segment_dict.keys():
        bolo = Bolo(bolo_name, config)
        for segment in bolo_segment_dict[bolo_name]:
            start_seg = time.time()
            prompter.prompt("Doing Bolo : %s Segment : %d Rank : %d" % (bolo_name, segment+1, rank))
            if "hitmap" in config.timestream_data_products:
                hitmap_local += bolo.simulate_timestream(segment)
            else:
                bolo.simulate_timestream(segment)
            stop_seg = time.time()
            prompter.prompt("Rank : " + str(rank) + " Time taken : " + str(stop_seg - start_seg) + ". Projected time : " + str((stop_seg - start_seg)*tot_seg))

    prompter.prompt("Done simulating")

    
    if "hitmap" in config.timestream_data_products:
        hitmap = np.zeros(hitmap_local.size)

        comm.Reduce(hitmap_local, hitmap, MPI.SUM, 0)
        
        scan_dir = os.path.join(config.general_data_dir, config.sim_tag, config.scan_tag)
        if rank == 0:
            hp.write_map(os.path.join(scan_dir, "hitmap_in.fits"), hitmap)



def make_data_dirs():
    sim_dir = os.path.join(config.general_data_dir, config.sim_tag)
    if not os.path.exists(sim_dir):
        try:
            os.makedirs(sim_dir)
        except OSError:
            pass
    scan_dir = os.path.join(sim_dir, config.scan_tag)
    if not os.path.exists(scan_dir):
        os.makedirs(scan_dir)
    config_dir = os.path.join(scan_dir, "config_files")
    if not os.path.exists(config_dir):
        os.makedirs(config_dir)
    #default_config_file = "/global/homes/b/banerji/simulation/timestream_simulation/config_files/default_config.py"
    #shutil.copy(default_config_file, config_dir)
    #current_config_file = os.path.join("/global/homes/b/banerji/simulation/timestream_simulation/config_files", config_file + '.py') 
    #shutil.copy(current_config_file, config_dir)


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
