#!/usr/bin/env python

import numpy as np
import os
import sys
import importlib
import time
from mpi4py import MPI
from simulation.lib.utilities.prompter import prompt

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

config_file = sys.argv[1]
config = importlib.import_module(config_file).config

num_segments = int(sys.argv[2])

sim_dir = os.path.join(config.general_output_dir, config.sim_tag)
scan_dir = os.path.join(sim_dir, config.scan_tag)

for bolo_pair in config.bolo_list:
    bolo_a_dir = os.path.join(scan_dir, bolo_pair + 'a')
    bolo_b_dir = os.path.join(scan_dir, bolo_pair + 'b')
    bolo_re_diff_dir = os.path.join(scan_dir, bolo_pair + 'a')
    gradient_co_dir = os.path.join(scan_dir, bolo_pair + 'a')
    gradient_cross_dir = os.path.join(scan_dir, bolo_pair + 'a')
    gradient_co_rescan_dir = os.path.join(scan_dir, bolo_pair + 'a')
    gradient_cross_rescan_dir = os.path.join(scan_dir, bolo_pair + 'a')

    local_segment_list = rank*num_segments/size + np.arange(num_segments/size)

    comm.Barrier()
    time.sleep(0.1*rank)
    prompt("Rank : {} doing segments : {}\n".format(rank, local_segment_list))
    comm.Barrier()

    local_num = 0
    local_den = 0

    for segment in local_segment_list:
        segment_name = str(segment+1).zfill(4)
        prompt("Rank : {} Doing segment : {}\n".format(rank, segment_name))
        signal_a = np.load(os.path.join(bolo_a_dir, segment_name, "signal.npy"))
        if np.sum(np.isnan(signal_a)):
            prompt("Rank : {}, Segment : {} in signal_a has nan\n".format(rank, segment_name))
        signal_b = np.load(os.path.join(bolo_b_dir, segment_name, "signal.npy"))
        if np.sum(np.isnan(signal_b)):
            prompt("Rank : {}, Segment : {} in signal_b has nan\n".format(rank, segment_name))
        signal_diff = 0.5*(signal_a - signal_b)
        signal_re_diff = np.load(os.path.join(bolo_re_diff_dir, segment_name, "resim_diff_signal.npy"))
        if np.sum(np.isnan(signal_re_diff)):
            prompt("Rank : {}, Segment : {} in signal_re_diff has nan\n".format(rank, segment_name))
        signal_TEMPLATE = np.empty((signal_a.size, 2))
        signal_re_TEMPLATE = np.empty((signal_a.size, 2))
        signal_TEMPLATE[:, 0] = np.load(os.path.join(gradient_co_dir, segment_name, "tm_grad_co.npy"))
        if np.sum(np.isnan(signal_TEMPLATE[:,0])):
            prompt("Rank : {}, Segment : {} in tm_grad_co has nan\n".format(rank, segment_name))
        signal_TEMPLATE[:, 1] = np.load(os.path.join(gradient_cross_dir, segment_name, "tm_grad_cross.npy"))
        if np.sum(np.isnan(signal_TEMPLATE[:,1])):
            prompt("Rank : {}, Segment : {} in tm_grad_cross has nan\n".format(rank, segment_name))
        signal_re_TEMPLATE[:, 0] = np.load(os.path.join(gradient_co_rescan_dir, segment_name, "resim_tm_grad_co.npy"))
        if np.sum(np.isnan(signal_re_TEMPLATE[:,0])):
            prompt("Rank : {}, Segment : {} in tm_resim_grad_co has nan\n".format(rank, segment_name))
        signal_re_TEMPLATE[:, 1] = np.load(os.path.join(gradient_cross_rescan_dir, segment_name, "resim_tm_grad_cross.npy"))
        if np.sum(np.isnan(signal_re_TEMPLATE[:,1])):
            prompt("Rank : {}, Segment : {} in tm_resim_grad_cross has nan\n".format(rank, segment_name))

        local_num += np.dot(signal_TEMPLATE.T, (signal_diff - signal_re_diff))
        local_den += np.dot(signal_TEMPLATE.T, (signal_TEMPLATE - signal_re_TEMPLATE))

    local_y = np.dot(np.linalg.inv(local_den), local_num)
    comm.Barrier()
    time.sleep(0.1*rank)

    prompt("Rank : {}\nLocal Numerator :\n{}\nLocal Denominator :\n{}\nLocal estimate of y : {}\n".format(rank, local_num, local_den, local_y))

    num = np.zeros(local_num.shape) 
    den = np.zeros(local_den.shape) 

    comm.Reduce(local_num, num, MPI.SUM, 0)
    comm.Reduce(local_den, den, MPI.SUM, 0)

    if rank == 0:
        np.save(os.path.join(sim_dir, bolo_pair + "_num_matrix.npy"), num)
        np.save(os.path.join(sim_dir, bolo_pair + "_den_matrix.npy"), den)
        y_global = np.dot(np.linalg.inv(den), num)
        prompt("Estimated y : {}\n".format(y_global)) 
        np.save(os.path.join(sim_dir, bolo_pair + "_estimated_y.npy"), y_global)
