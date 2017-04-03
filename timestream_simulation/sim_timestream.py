#!/usr/bin/env python

import numpy as np
import healpy as hp
import sys
import os
import shutil
import importlib
import time
import pickle
from memory_profiler import profile
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list
from simulation.lib.utilities.prompter import prompt
from simulation.timestream_simulation.bolo import Bolo


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* The master section which distributes the data and collects it
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*


def run_simulation():
    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)
    tot_seg = 0
    for keys in bolo_segment_dict.keys():
        tot_seg += len(bolo_segment_dict[keys])

    if run_type == 'run_mpi':
        comm.Barrier()
        time.sleep(0.1*rank)
        prompt("Rank : {} Local bolo segment list :\n{}\n# of segments : {}\n".format(rank, bolo_segment_dict, tot_seg), sys.stdout)
        comm.Barrier()
    else:
        prompt("Rank : {} Local bolo segment list :\n{}\n# of segments : {}\n".format(rank, bolo_segment_dict, tot_seg), sys.stdout)

    comm.Barrier()

    count = 0
    time_taken = 0

    for bolo_name in bolo_segment_dict.keys():
        bolo = Bolo(bolo_name, config)
        for segment in bolo_segment_dict[bolo_name]:
            count += 1
            start_seg = time.time()
            prompt("Doing Bolo : {} Segment : {} Rank : {}\n".format(bolo_name, segment+1, rank), sys.stdout)
            if config.sim_type == "signal":
                bolo.simulate_timestream_signal(segment)
            else:
                bolo.simulate_timestream_template(segment)
            stop_seg = time.time()
            time_taken += stop_seg - start_seg
            prompt("Rank : {}, Time taken : {}. Total time take : {}, Projected time : {}, Finished {} of {}\n".format(rank, stop_seg - start_seg, time_taken, time_taken*tot_seg/count, count, tot_seg))

    prompt("Rank : {} Done simulating\n".format(rank), sys.stdout)


def make_data_dirs():
    sim_dir = os.path.join(config.general_output_dir, config.sim_tag)
    if not os.path.exists(sim_dir):
        try:
            os.makedirs(sim_dir)
        except OSError:
            pass

    scan_dir = os.path.join(sim_dir, config.scan_tag)
    if not os.path.exists(scan_dir):
        try:
            os.makedirs(scan_dir)
        except OSError:
            pass

    config_dir = os.path.join(sim_dir, "config_files")
    if not os.path.exists(config_dir):
        try:
            os.makedirs(config_dir)
        except OSError:
            pass
    this_config_dir = os.path.join(config_dir, time_stamp)
    try:
        os.makedirs(this_config_dir)
    except OSError:
        pass
    with open(os.path.join(this_config_dir, "config_file.pkl"), "w") as outfile:
        pickle.dump(config, outfile)
    bolo_config = importlib.import_module(config.bolo_config_file).bolo_config
    with open(os.path.join(this_config_dir, "bolo_config_file.pkl"), "w") as outfile:
        pickle.dump(bolo_config, outfile)

    code_dir = os.path.join(sim_dir, "code_files")
    if not os.path.exists(code_dir):
        try:
            os.makedirs(code_dir)
        except OSError:
            pass
    this_code_dir = os.path.join(code_dir, time_stamp)
    try:
        os.makedirs(this_code_dir)
    except OSError:
        pass
    shutil.copyfile(os.path.join(config.base_dir, "timestream_simulation", "sim_timestream.py"), os.path.join(this_code_dir, "sim_timestream.py"))
    shutil.copyfile(os.path.join(config.base_dir, "timestream_simulation", "bolo.py"), os.path.join(this_code_dir, "bolo.py"))
    shutil.copyfile(os.path.join(config.base_dir, "timestream_simulation", "beam_kernel.py"), os.path.join(this_code_dir, "beam_kernel.py"))


def start_message():
    display_string = "\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"
    display_string += "#* BEGINNING SIMULATION\n"
    display_string += "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"
    display_string += "TIME STAMP : {}\n".format(time_stamp)
    display_string += "RUN TYPE : {}\n".format(run_type)
    display_string += "# OF PROCESSES : {}\n".format(size)
    display_string += "DETECTOR LIST : {}\n".format(config.bolo_list)
    display_string += "# OF DETECTORS : {}\n".format(len(config.bolo_list))
    display_string += "SEGMENT LIST : {}\n".format(config.segment_list)
    display_string += "# OF SEGMENTS : {}\n".format(len(config.segment_list))
    display_string += "SIMULATION TYPE : {}\n".format(config.sim_type)
    display_string += "BEAM TYPE : {}\n".format(config.beam_type)
    if config.sim_type == "template":
        if config.template_type == "tm_gradient":
            display_string += "GRADIENT TYPES : {}\n".format(config.tm_gradient_type)
    display_string += "WRITE FIELD : {}\n".format(config.timestream_data_products)
    display_string += "NOTES : {}\n".format(config.notes)
    display_string += "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n"
    prompt(display_string)

    for bolo_name in config.bolo_list:
        bolo = Bolo(bolo_name, config)
        bolo.display_params()
        bolo.beam.display_beam_settings()

def run_check():
    if config.sim_type == "gradient" and config.sim_pol_type != 'T':
        if rank == 0:
            prompt("When computing gradients, the simulation polarisation is set to T only. Changing", color="red")
        config.sim_pol_type = 'T'

    if config.sim_pol_type == "noise_only" and config.noise_type == "none":
        if rank == 0:
            prompt("You need to provide the type of noise when doing a noise only simulation. Otherwise the output will be just zero for all array elements.", color="red")
        sys.exit()

    if config.noise_type == "1_over_f":
        if rank == 0:
            prompt("One over f noise is not developed yet. Please be patient.", color="red")
        sys.exit()


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#* Main function definition. This is where the code begins when executed
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

if __name__=="__main__":
    config_file = sys.argv[1]
    run_type = sys.argv[2]

    config = importlib.import_module(config_file).config

    time_stamp = get_time_stamp()

    if run_type == 'run_mpi':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        comm.bcast(time_stamp, root=0)
        if rank == 0:
            make_data_dirs()
            start_message()
        comm.Barrier()

    if run_type == 'run_serial':
        size = 1
        rank = 0
        make_data_dirs()
        start_message()
    
    run_check()

    run_simulation()
