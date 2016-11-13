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

data_dir = config.general_data_dir
sim_dir = os.path.join(data_dir, config.sim_tag)
bolo_a_dir = os.path.join(sim_dir, "first_scan", "bolo_0001a") 
bolo_b_dir = os.path.join(sim_dir, "first_scan", "bolo_0001b") 
bolo_TEMPLATE_dir = os.path.join(sim_dir, "first_scan", "bolo_TEMPLATE")

num_segments = int(sys.argv[2])

local_segment_list = rank*num_segments/size + np.arange(num_segments/size)
#print "Rank :", rank, "Local segment list :", local_segment_list

local_num = 0
local_den = 0

#sys.exit()

for segment in local_segment_list:
    segment_name = str(segment+1).zfill(4)
    print "Rank :", rank, "Doing segment :", segment_name
    signal_a = np.load(os.path.join(bolo_a_dir, segment_name, "signal.npy"))
    signal_b = np.load(os.path.join(bolo_b_dir, segment_name, "signal.npy"))
    signal_diff = 0.5*(signal_a - signal_b)
    signal_TEMPLATE = np.load(os.path.join(bolo_TEMPLATE_dir, segment_name, "signal.npy"))
    local_num += np.sum(signal_TEMPLATE*signal_diff)
    local_den += np.sum(signal_TEMPLATE**2)

comm.Barrier()

print "Rank :", rank, "Local estimate of y :", local_num/local_den

local_num = np.array(local_num)
local_den = np.array(local_den)

num = np.zeros(1) 
den = np.zeros(1) 

comm.Reduce(local_num, num, MPI.SUM, 0)
comm.Reduce(local_den, den, MPI.SUM, 0)

if rank == 0:
    print "Estimated y :", num[0]/den[0]
    print "Estimated delta alpha :", 4*num[0]/den[0]
    np.save("estimated_y", num[0]/den[0])
