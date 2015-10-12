import numpy as np
import os
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.global_settings import global_paths, global_scanning

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Parameters 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

settings.nside = global_scanning.nside
settings.load_pointing = False
settings.generate_pointing = True
settings.do_pol = True
settings.do_beam_kernel = False 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

tag = "test"
base_folder = global_paths.base_folder 
maps_folder = global_paths.maps_folder
settings.input_map = os.path.join(maps_folder, "test_map_2048_0.fits")

settings.output_folder = global_paths.output_folder 

settings.write_signal = False 
settings.display_scanned_map = True 
settings.write_scanned_map = True
settings.pipe_with_map_maker = False 
