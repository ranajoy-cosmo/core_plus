import numpy as np
from simulation.lib.utilities.generic_class import Generic
from simulation.params.custom_params import global_paths, global_scanning
from simulation.lib.utilities.time_util import get_time_stamp
import os

map_making_params = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Resolution params
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
map_making_params.nside_out = global_scanning.nside_out

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Params 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
map_making_params.time_stamp = get_time_stamp()
map_making_params.scanning_time_stamp = "" 
map_making_params.bolo_names = ["0001", "0002", "0003", "0004"]
map_making_params.t_flight = 4*24*60*60.0
map_making_params.t_segment = 12*60*60.0
map_making_params.global_output_dir = global_paths.output_dir
map_making_params.display_map = False 
map_making_params.write_map = True
