#!/usr/bin/env python

import numpy as np
import os
import sys
import importlib
from mpi4py import MPI
from simulation.lib.prompter import prompt

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

config_file = sys.argv[1]
config = importlib.import_module(config_file).config

num_segments = int(sys.argv[2])

data_dir = config.general_data_dir
sim_dir = os.path.join(data_dir, config.sim_tag)

for bolo_pair in config.bolo_list:
    bolo_a_dir = os.path.join(sim_dir, "scan", bolo_pair + "a") 
    bolo_b_dir = os.path.join(sim_dir, "scan", bolo_pair + "b") 
    bolo_re_diff_dir = os.path.join(sim_dir, "scan", bolo_pair + "_diff_re")
    TEMPLATE_dir_list = []
    TEMPLATE_re_dir_list = []
    for TEMPLATE_dir in config.bolo_TEMPLATE_dir_dict[bolo_pair]:
        TEMPLATE_dir_list.append(os.path.join(sim_dir, "scan", bolo_pair + '_' + TEMPLATE_dir))
        TEMPLATE_re_dir_list.append(os.path.join(sim_dir, "scan", bolo_pair + '_' + TEMPLATE_dir + 're'))

    local_segment_list = rank*num_segments/size + np.arange(num_segments/size)

    local_num = 0
    local_den = 0

    for segment in local_segment_list:
        segment_name = str(segment+1).zfill(4)
        print "Rank : {} Doing segment : {}\n".format(rank, segment_name))
        signal_a = np.load(os.path.join(bolo_a_dir, segment_name, "signal.npy"))
        signal_b = np.load(os.path.join(bolo_b_dir, segment_name, "signal.npy"))
        signal_diff = 0.5*(signal_a - signal_b)
        signal_re_diff = np.load(os.path.join(bolo_re_diff_dir, segment_name, "signal.npy"))
        signal_TEMPLATE = np.empty((signal_a.size, len(TEMPLATE_dir_list)))
        signal_re_TEMPLATE = np.empty((signal_a.size, len(TEMPLATE_dir_list)))
        for i in range(len(TEMPLATE_dir_list)):
            signal_TEMPLATE[..., i] = np.load(os.path.join(TEMPLATE_dir_list[i], segment_name, "signal.npy"))
            signal_re_TEMPLATE[..., i] = np.load(os.path.join(TEMPLATE_re_dir_list[i], segment_name, "signal.npy"))

        local_num += np.dot(signal_TEMPLATE.T, (signal_diff - signal_re_diff))
        local_den += np.dot(signal_TEMPLATE.T, (signal_TEMPLATE - signal_re_TEMPLATE))

    comm.Barrier()

    local_y = np.dot(np.linalg.inv(local_den), local_num)
    print "Rank : {} Local estimate of y : {}\n".format(rank, local_y))

    num = np.zeros(local_num.shape) 
    den = np.zeros(local_den.shape) 

    comm.Reduce(local_num, num, MPI.SUM, 0)
    comm.Reduce(local_den, den, MPI.SUM, 0)

    if rank == 0:
        np.save(os.path.join(sim_dir, bolo_pair + "_num_matrix.npy", num))
        np.save(os.path.join(sim_dir, bolo_pair + "_den_matrix.npy", den))
        y_global = np.dot(np.linalg.inv(den), num)
        prompt("Estimated y : {}\n".format(y_global)) 
        np.save(os.path.join(config.general_data_dir, config.sim_tag, bolo_pair + "_estimated_y.npy"), y_global)
