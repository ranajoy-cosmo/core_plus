#! /usr/bin/env python

import numpy as np
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.custom_settings import global_paths, global_scanning
import os

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Map resolution settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.lmax = 3950
settings.nside = global_scanning.nside_in
settings.nside_deg = global_scanning.nside_out

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam and Polarisation Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.fwhm_arcmin = 8.0                                          #arc-minutes
settings.fwhm = np.deg2rad(settings.fwhm_arcmin/60.0)               #radians
settings.do_unbeamed_map = True
settings.do_beamed_map = True 
settings.do_degraded_map = True 
settings.do_pixel_window = True
settings.do_pol = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.tag = 'sky'
spectra_file_name = "test_4000.npy"
settings.spectra_file = os.path.join(global_paths.spectra_dir, spectra_file_name) 

map_file_name = settings.tag + '_map_'
settings.map_file = os.path.join(global_paths.maps_dir, map_file_name) 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# I/O Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.write_map = True
settings.return_map = False
settings.display_maps = False
