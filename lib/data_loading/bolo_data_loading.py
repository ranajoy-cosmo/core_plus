#!/usr/bin/env python 

import numpy as np
import healpy as hp
import os
import sys

def check_data_request_validity(scan_params):
    t_flight = scan_params.t_flight
    t_segment = scan_params.

def create_output_sim_dirs(time_stamp, scan_params):
    out_dir = os.path.join(scan_params.global_output_dir, "scanning", time_stamp)
