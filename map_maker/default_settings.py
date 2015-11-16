import numpy as np
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.custom_settings import global_paths, global_scanning
import simulation.bolo.custom_settings as scan_settings
import os

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Resolution params
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.nside_in = global_scanning.nside_in
settings.nside_out = global_scanning.nside_out

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Params 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
base_folder = global_paths.base_folder
settings.output_folder = global_paths.output_folder
settings.display_map = False 
settings.write_map = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Params 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
scan_params = scan_settings.settings
scan_params.write_signal = True
scan_params.pipe_with_map_maker = False
scan_params.display_scanned_map = False
