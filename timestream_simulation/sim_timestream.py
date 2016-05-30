#!/usr/bin/rnv python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
import os
import simulation.lib.data_management.data_utilities as data_utils
from simulation.lib.utilities.time_util import get_time_stamp

class Bolo:

    def __init__(self, bolo_name):
        self.bolo_params = importlib.import_module("simulation.timestream_simulation.bolo_params." + bolo_name).bolo_params
        self.beam_kernel, self.del_beta = 
        self.sky_map = self.load_sky_map()
        self.pointing = Pointing(self.bolo_params, scan_params)



def run_mpi():
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if scan_params.tag == None:
        tag = [get_time_stamp()]
        comm.bcast(tag, root=0)
        tag = tag[0]
    else:
        tag = scan_params.tag

    out_dir = os.path.join(scan_params.global_output_dir, "scanning", tag)

    calculate_params()

    if rank == 0:
        display_params()
        data_utils.make_output_dirs(out_dir, scan_params.bolo_list, scan_params.segment_list)
        data_utils.copy_source_to_output(out_dir, scan_params.bolo_list)

    comm.Barrier()

    bolo_segment_dict = bolo_data_loading.get_local_bolo_segment_list(rank, size, scan_params)

    for bolo_name in bolo_segment_dict.keys():
        bolo = Bolo(bolo_name)
        for segment in bolo_segment_dict[bolo_name]:
            bolo.generate_timestream(segment)


if __name__=="__main__":
    time_start = time.time()
    from custom_params import scan_params, beam_params
    action = sys.argv[1]

    if action=='run_mpi':
        from mpi4py import MPI
        run_mpi()

    time_stop = time.time()
    print "Total time taken :", (time_stop - time_start), "s"
