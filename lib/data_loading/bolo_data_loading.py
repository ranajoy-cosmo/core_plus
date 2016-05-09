#!/usr/bin/env python 

import numpy as np

def get_local_bolo_segment_list(rank, num_processes, params):
    num_bolos = len(params.bolo_names)
    num_segments_per_bolo = params.segment_list.size
    num_total_segments = num_bolos*num_segments_per_bolo

    if num_total_segments % num_processes != 0:
        num_segments_per_process = num_total_segments/num_processes + 1
    else:
        num_segments_per_process = num_total_segments/num_processes
    
    bolo_list = []
    for bolo_name in params.bolo_names:
        bolo_list.extend([bolo_name]*num_segments_per_bolo)

    segment_list = np.arange(num_total_segments) % num_segments_per_bolo

    local_bolo_list = bolo_list[rank*num_segments_per_process : (rank + 1)*num_segments_per_process]
    local_segment_list = segment_list[rank*num_segments_per_process : (rank + 1)*num_segments_per_process]

    return local_bolo_list, local_segment_list

