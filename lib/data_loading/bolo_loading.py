#!/usr/bin/env python 

import numpy as np
import healpy as hp
import os, sys
from simulation.settings import custom_settings

def create_output_sim_dirs(data_params):
    out_dir = os.path.join
