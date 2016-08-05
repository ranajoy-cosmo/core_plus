import numpy as np
from simulation.lib.utilities.generic_class import Generic
from simulation.params.custom_params import global_paths, global_scanning
from simulation.lib.utilities.time_util import get_time_stamp
import os


config = Generic()


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Run
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# config.run_type = ["serial"/"mpi"/"display_parameters"]
config.run_type = "mpi"
#config.simulate_ts = True/False. True means it will simulate the entire timestream signal and pointing from scratch
config.simulate_ts = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Resolution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.nside_out : nside of the output map
config.nside_out = 1024 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.global_output_dir = global_paths.output_dir
config.tag = get_time_stamp()
config.write_output_map = True
config.write_output_hitmap = True
config.write_inverse_covariance_matrix = True
config.write_covariance_matrix = True
config.write_inverse_pol_covariance_matrix = True
config.write_pol_covariance_matrix = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#This section is only valid when config.simulate_ts is False

config.scanning_tag = "" 
config.bolo_list = ["bolo_0001a"]
config.segment_list = np.arange(1)


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulation configurations 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#This section is only valid when config.simulate_ts is True 

import simulation.timestream_simulation.default_config.config as sim_config
config.__dict__.update(sim_config.__dict__)
