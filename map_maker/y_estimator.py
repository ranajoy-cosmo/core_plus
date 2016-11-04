#!/usr/bin/env python

import numpy as np
import os
import sys
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

sim_dir = "/home/banerji/simulation/output/bandpass_test_gal_6"
bolo_a_dir = os.path.join(sim_dir, "first_scan", "bolo_0001a") 
bolo_b_dir = os.path.join(sim_dir, "first_scan", "bolo_0001b") 
bolo_TEMPLATE_dir = os.path.join(sim_dir, "first_scan", "bolo_TEMPLATE")
bolo_diff_rescan_dir = os.path.join(sim_dir, "rescan", "bolo_0001a") 
bolo_TEMPLATE_rescan_dir = os.path.join(sim_dir, "rescan", "bolo_TEMPLATE") 

#print os.path.exists(bolo_a_dir)
#print os.path.exists(bolo_TEMPLATE_dir)
#print os.path.exists(bolo_b_dir)
#print os.path.exists(bolo_TEMPLATE_dir)
#print os.path.exists(bolo_diff_rescan_dir)
#print os.path.exists(bolo_TEMPLATE_rescan_dir)

num_segments = int(sys.argv[1])

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
    signal_diff_rescan = np.load(os.path.join(bolo_diff_rescan_dir, segment_name, "signal.npy"))
    signal_TEMPLATE = np.load(os.path.join(bolo_TEMPLATE_dir, segment_name, "signal.npy"))
    signal_TEMPLATE_rescan = np.load(os.path.join(bolo_TEMPLATE_rescan_dir, segment_name, "signal.npy"))
    local_num += np.sum(signal_TEMPLATE*(signal_diff - signal_diff_rescan))
    local_den += np.sum(signal_TEMPLATE*(signal_TEMPLATE - signal_TEMPLATE_rescan))

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
