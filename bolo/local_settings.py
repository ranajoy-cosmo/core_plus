import numpy as np
import os
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.global_settings import global_paths, global_scanning
import simulation.pointing.local_settings as pointing_settings

scan_params = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Parameters 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

scan_params.nside = global_scanning.nside
scan_params.load_pointing = False
scan_params.generate_pointing = True
scan_params.do_pol = True
scan_params.do_beam_kernel = False 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Pointing Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

pointing_params = pointing_settings.settings
pointing_params = 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

rw_params = Generic()
rw_params.tag = "test"
rw_params.base_folder = global_paths.base_folder 
rw_params.maps_folder = global_paths.maps_folder
rw_params.input_map = os.path.join(maps_folder, "test_map_2048_0.fits")

rw_params.output_folder = global_paths.output_folder 

rw_params.write_signal = False 
rw_params.display_scanned_map = True 
rw_params.settings.write_scanned_map = True
rw_params.settings.pipe_with_map_maker = False 
