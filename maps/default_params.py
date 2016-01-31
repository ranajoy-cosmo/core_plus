#! /usr/bin/env python

import numpy as np
from simulation.lib.utilities.generic_class import Generic
from simulation.params.custom_params import global_paths, global_scanning
import os

map_params = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Map resolution settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
map_params.lmax = 3950
map_params.nside = global_scanning.nside_in
map_params.nside_deg = global_scanning.nside_out

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam and Polarisation Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
map_params.fwhm = 8.0                                          #arc-minutes
map_params.do_unbeamed_map = True
map_params.do_beamed_map = True 
map_params.do_degraded_map = True 
map_params.do_pixel_window = True
map_params.do_pol = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
map_params.tag = 'sky'
spectra_file_name = "test_4000.npy"
map_params.spectra_file = os.path.join(global_paths.spectra_dir, spectra_file_name) 

map_file_name = map_params.tag + '_map_'
map_params.map_file = os.path.join(global_paths.maps_dir, map_file_name) 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# I/O Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
map_params.write_map = True
map_params.return_map = False
map_params.display_maps = False
