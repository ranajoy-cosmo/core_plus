import numpy as np
import os
import sys
import shutil
import simulation.lib.utilities.prompter as prompter

def get_local_bolo_segment_list(rank, num_processes, bolo_list, segment_list):
    num_bolos = len(bolo_list)
    num_segments_per_bolo = len(segment_list)
    num_total_segments = num_bolos*num_segments_per_bolo

    num_segments_per_process = num_total_segments/num_processes

    bolo_list_deg = []
    for bolo_name in bolo_list:
        bolo_list_deg.extend([bolo_name]*num_segments_per_bolo)

    segment_list_deg = list(np.arange(num_total_segments) % num_segments_per_bolo)

    local_bolo_list = bolo_list_deg[rank*num_segments_per_process : (rank + 1)*num_segments_per_process]
    local_segment_list = segment_list_deg[rank*num_segments_per_process : (rank + 1)*num_segments_per_process]

    if num_total_segments % num_processes != 0:
        num_extra_segments = num_total_segments%num_processes
        extra_bolo_list = bolo_list_deg[-num_extra_segments:]
        extra_segment_list = segment_list_deg[-num_extra_segments:]
        if rank<num_extra_segments:
            local_bolo_list.extend([extra_bolo_list[rank]])
            local_segment_list.extend([extra_segment_list[rank]])

    bolo_segment_dict = dict.fromkeys(set(local_bolo_list), None)

    for bolo in bolo_segment_dict.keys():
        bolo_segment_dict[bolo] = [local_segment_list[i] for i in xrange(len(local_segment_list)) if local_bolo_list[i] == bolo]

    return bolo_segment_dict


def get_bolo_pair_segment_list(rank, num_processes, segment_list):
    num_segments = len(segment_list)

    if num_segments % num_processes != 0:
        num_segments_per_process = num_segments/num_processes + 1
    else:
        num_segments_per_process = num_segments/num_processes

    pair_segment_list = segment_list[rank*num_segments_per_process : (rank+1)*num_segments_per_process]

    return pair_segment_list
