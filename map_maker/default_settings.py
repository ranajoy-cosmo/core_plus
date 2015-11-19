import numpy as np
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.custom_settings import global_paths, global_scanning
import simulation.bolo.custom_settings as scan_settings
from simulation.lib.utilities.time_util import get_time_stamp
import os

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Resolution params
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.nside_out = global_scanning.nside_out

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Params 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.time_stamp = get_time_stamp()
settings.scanning_time_stamp = "19_36_33__19_11_2015"
settings.bolo_names = ["0001", "0002"]
settings.t_flight = 10*60*60.0
settings.t_segment = 5*60*60.0
settings.global_output_dir = global_paths.output_dir
settings.display_map = False 
settings.write_map = True
