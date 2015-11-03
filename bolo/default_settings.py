import numpy as np
import os
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.global_settings import global_paths, global_scanning
import simulation.pointing.default_settings as pointing_settings

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Parameters 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

settings.bolo_names = ['0001']
settings.nside_in = global_scanning.nside_in
settings.nside_out = global_scanning.nside_out
settings.load_pointing = False
settings.do_pol = False

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Pointing Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

pointing_params = pointing_settings.settings
pointing_params.write_pointing = True
pointing_params.return_pointing = True
pointing_params.display_params = False
pointing_params.do_pol = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

settings.base_folder = global_paths.base_folder 
settings.input_map = os.path.join(global_paths.maps_folder, "test_map_1024_0.fits")
settings.bolo_param_folder = os.path.join(global_paths.base_folder, "bolo", "bolo_params")
settings.output_folder = global_paths.output_folder 

settings.write_signal = True 
settings.display_scanned_map = True 
settings.write_scanned_map = True
settings.write_hitmap = True
settings.pipe_with_map_maker = False 
