#!/usr/bin/env python

import numpy as np
import os
import sys
import importlib
from mpi4py import MPI

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
    bolo_TEMPLATE_TD_dir = os.path.join(sim_dir, "scan", bolo_pair + "_TEMPLATE_TD")
    bolo_TEMPLATE_SYNC_dir = os.path.join(sim_dir, "scan", bolo_pair + "_TEMPLATE_SYNC")
    bolo_re_TEMPLATE_TD_dir = os.path.join(sim_dir, "scan", bolo_pair + "_TEMPLATE_TD_re")
    bolo_re_TEMPLATE_SYNC_dir = os.path.join(sim_dir, "scan", bolo_pair + "_TEMPLATE_SYNC_re")

    local_segment_list = rank*num_segments/size + np.arange(num_segments/size)
    #print "Rank :", rank, "Local segment list :", local_segment_list

    local_num = 0
    local_den = 0

    for segment in local_segment_list:
        segment_name = str(segment+1).zfill(4)
        print "Rank :", rank, "Doing segment :", segment_name
        signal_a = np.load(os.path.join(bolo_a_dir, segment_name, "signal.npy"))
        signal_b = np.load(os.path.join(bolo_b_dir, segment_name, "signal.npy"))
        signal_diff = 0.5*(signal_a - signal_b)
        signal_re_diff = np.load(os.path.join(bolo_re_diff_dir, segment_name, "signal.npy"))
        signal_TEMPLATE = np.empty((signal_a.size, 2))
        signal_TEMPLATE[..., 0] = np.load(os.path.join(bolo_TEMPLATE_TD_dir, segment_name, "signal.npy"))
        signal_TEMPLATE[..., 1] = np.load(os.path.join(bolo_TEMPLATE_SYNC_dir, segment_name, "signal.npy"))
        signal_re_TEMPLATE = np.empty((signal_a.size, 2))
        signal_re_TEMPLATE[..., 0] = np.load(os.path.join(bolo_re_TEMPLATE_TD_dir, segment_name, "signal.npy"))
        signal_re_TEMPLATE[..., 1] = np.load(os.path.join(bolo_re_TEMPLATE_SYNC_dir, segment_name, "signal.npy"))

        local_num += np.dot(signal_TEMPLATE.T, (signal_diff - signal_re_diff))
        local_den += np.dot(signal_TEMPLATE.T, (signal_TEMPLATE - signal_re_TEMPLATE))

    comm.Barrier()

    local_y = np.dot(np.linalg.inv(local_den), local_num)
    print "Rank :", rank, "Local estimate of y :", local_y

    num = np.zeros(local_num.shape) 
    den = np.zeros(local_den.shape) 

    comm.Reduce(local_num, num, MPI.SUM, 0)
    comm.Reduce(local_den, den, MPI.SUM, 0)

    if rank == 0:
        y_global = np.dot(np.linalg.inv(den), num)
        print "Estimated y :", y_global 
        np.save(os.path.join(config.general_data_dir, config.sim_tag, bolo_pair + "_estimated_y.npy"), y_global)
        print "Writing at :", os.path.join(config.general_data_dir, config.sim_tag, bolo_pair + "_estimated_y.npy")
