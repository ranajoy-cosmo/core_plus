#!/usr/bin/env python 

import numpy as np
import os
import shutil

def get_local_bolo_segment_list(rank, num_processes, bolo_list, segment_list):
    num_bolos = len(bolo_list)
    num_segments_per_bolo = segment_list.size
    num_total_segments = num_bolos*num_segments_per_bolo

    if num_total_segments % num_processes != 0:
        num_segments_per_process = num_total_segments/num_processes + 1
    else:
        num_segments_per_process = num_total_segments/num_processes

    bolo_list_deg = []
    for bolo_name in bolo_list:
        bolo_list_deg.extend([bolo_name]*num_segments_per_bolo)

    segment_list_deg = np.arange(num_total_segments) % num_segments_per_bolo

    local_bolo_list = bolo_list_deg[rank*num_segments_per_process : (rank + 1)*num_segments_per_process]
    local_segment_list = segment_list_deg[rank*num_segments_per_process : (rank + 1)*num_segments_per_process]

    bolo_segment_dict = dict.fromkeys(set(local_bolo_list), None)

    for bolo in bolo_segment_dict.keys():
        bolo_segment_dict[bolo] = [local_segment_list[i] for i in xrange(len(local_segment_list)) if local_bolo_list[i] == bolo]

    return bolo_segment_dict

def create_output_directories(output_dir, bolo_list, segment_list):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    for bolo in bolo_list:
        if not os.path.isdir(os.path.join(output_dir, bolo)):
            os.makedirs(os.path.join(output_dir, bolo))

        for segment in segment_list:
            segment_name = str(segment+1).zfill(4)
            if not os.path.isdir(os.path.join(output_dir, bolo, segment_name)):
                os.makedirs(os.path.join(output_dir, bolo, segment_name))

def copy_source_to_output(output_dir, bolo_list):
    param_dir = os.path.join(output_dir, "params")
    if not os.path.isdir(param_dir):
        os.makedirs(param_dir)

    shutil.copy("default_params.py", param_dir)
    shutil.copy("custom_params.py", param_dir)
    shutil.copy("sim_timestream_pol.py", output_dir)

    bolo_param_dir = os.path.join(param_dir, "bolo_params")
    if not os.path.isdir(bolo_param_dir):
        os.makedirs(bolo_param_dir)
    for bolo in bolo_list:
        shutil.copy("bolo_params/"+bolo+".py", bolo_param_dir)
