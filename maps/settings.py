#! /usr/bin/env python

import numpy as np
from simulation.lib.utilities.generic_class import Generic

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Map resolution settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.lmax = 4000
settings.nside = 2048 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam and Polarisation Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.fwhm_arcmin = 8.0
settings.fwhm = np.deg2rad(settings.fwhm_arcmin/60.0)
settings.do_unbeamed_map = True
settings.do_beamed_map = True 
settings.do_pol = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.tag = "test_map"
root_folder = "/Users/banerji/CORE+/simulation/"

spectrum_input_folder = root_folder + "maps_and_spectra/spectra/"
spectra_file_name = "test_4000.npy"
settings.spectra_file = spectrum_input_folder + spectra_file_name  

output_map_folder = root_folder + "maps_and_spectra/maps/"
map_file_name = settings.tag + '_' + str(settings.nside)
settings.map_file = output_map_folder + map_file_name

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# I/O Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.write_map = True
settings.return_map = False
settings.display_maps = False
