#! /usr/bin/env python
import numpy as np
from simulation.lib.utilities.generic_class import Generic

settings = Generic()

settings.tag = "test"

settings.lmax = 3000

settings.nside = 1024 

fwhm_arcmin = 8.0

settings.fwhm = np.deg2rad(fwhm_arcmin/60.0)

root_folder = "/Users/banerji/CORE+/simulation/"

spectrum_input_folder = root_folder + "maps_and_spectra/spectra/"
spectrum_file_name = "test_3000.npy"
settings.spectra_file = spectrum_input_folder + spectrum_file_name  

output_map_folder = root_folder + "maps_and_spectra/maps/"
map_file_name = settings.tag + '_' + str(settings.nside)
settings.map_file = output_map_folder + map_file_name

settings.write_map = True
settings.return_map = False
settings.display_maps = True
settings.do_pol = True

settings.do_beam_smoothing = False 
