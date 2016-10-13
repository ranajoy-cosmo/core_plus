#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import shutil
import time
import importlib
#from memory_profiler import profile
from simulation.lib.data_management.data_utilities import get_local_bolo_segment_list 
from simulation.timestream_simulation.bolo import Bolo
import simulation.lib.utilities.prompter as prompter
from simulation.lib.utilities.time_util import get_time_stamp 

def run_mpi():
    start = time.time()

    bolo_segment_dict = get_local_bolo_segment_list(rank, size, config.bolo_list, config.segment_list)

    print "Rank :", rank, ", Segments :", bolo_segment_dict
    comm.Barrier()

    num_bolos = len(config.bolo_list)
    num_local = np.zeros(num_bolos)
    den_local = np.zeros(num_bolos)

    for bolo_name in bolo_segment_dict.keys():
        bolo_name_a = bolo_name + 'a'
        bolo_name_b = bolo_name + 'b'
        for segment in bolo_segment_dict[bolo_name]:
            start_segment = time.time()
            #Bolo 145
            bolo_a = Bolo(bolo_name_a, config)
            bolo_b = Bolo(bolo_name_b, config)
            prompter.prompt("Rank : %d doing Bolo : %s and segment : %d" % (rank, bolo_name, segment))
            signal_b, v, pol_ang = bolo_b.read_timestream(segment)
            signal_a, v, pol_ang = bolo_a.read_timestream(segment)
            signal_diff = 0.5*(signal_a - signal_b)
            #Bolo 350
            bolo_350 = Bolo("bolo_00350", config)
            signal_350, v, pol_ang = bolo_350.read_timestream(segment) 
            bin_data(num_local, den_local, np.sum(signal_diff*signal_350), np.sum(signal_350*signal_350), bolo_name)
            stop_segment = time.time()
            time_per_seg = stop_segment - start_segment
            prompter.prompt("Time taken : %d" %(time_per_seg))


    num = np.zeros(num_bolos)
    den = np.zeros(num_bolos)

    comm.Reduce(num_local, num, MPI.SUM, 0)
    comm.Reduce(den_local, den, MPI.SUM, 0)

    if rank == 0:
        recon_dir = os.path.join(config.general_data_dir, config.sim_tag, config.scan_tag)
        alpha = dict(zip(config.bolo_list, num/den))
        print alpha
        np.save(os.path.join(recon_dir, "alpha"), alpha)
        
        stop = time.time()

        prompter.prompt("Total time taken : %d" % (stop - start))

def run_serial():

    num_bolos = len(config.bolo_list)
    num = np.zeros(num_bolos)
    den = np.zeros(num_bolos)

    for bolo_name in config.bolo_list:
        bolo_name_a = bolo_name + 'a'
        bolo_name_b = bolo_name + 'b'
        for segment in config.segment_list: 
            #Bolo 145
            bolo_a = Bolo(bolo_name_a, config)
            bolo_b = Bolo(bolo_name_b, config)
            prompter.prompt("Doing Bolo : %s and segment : %d" % (bolo_name, segment))
            signal_b, v, pol_ang = bolo_b.read_timestream(segment)
            signal_a, v, pol_ang = bolo_a.read_timestream(segment)
            signal_diff = signal_a - signal_b
            #Bolo 350
            prompter.prompt("Doing Bolo : %s and segment : %d" % ("bolo_350", segment))
            bolo_350 = Bolo("bolo_00350", config)
            signal_350, v, pol_ang = bolo_350.read_timestream(segment) 
            bin_data(num, den, np.sum(signal_diff*signal_350), np.sum(signal_350*signal_350), bolo_name)

    recon_dir = os.path.join(config.general_data_dir, config.sim_tag, config.scan_tag)
    alpha = dict(zip(config.bolo_list, num/den))
    print alpha
    np.save(os.path.join(recon_dir, "alpha"), alpha)


def bin_data(num, den, num_temp, den_temp, bolo_name):
    num[bolo_dict[bolo_name]] += num_temp
    den[bolo_dict[bolo_name]] += den_temp


if __name__=="__main__":
    config_file = sys.argv[1]
    run_type = sys.argv[2]

    config = importlib.import_module("simulation.map_maker.config_files." + config_file).config

    num_bolos = len(config.bolo_list)
    num_bin = np.zeros(num_bolos)
    den_bin = np.zeros(num_bolos)
    bolo_dict = dict(zip(config.bolo_list, range(num_bolos)))

    if run_type=='run_mpi':
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        run_mpi()

    if run_type=='run_serial':
        run_serial()
