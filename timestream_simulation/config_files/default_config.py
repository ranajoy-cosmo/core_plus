import healpy as hp
import os
import numpy as np
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.lib.utilities.generic_class import Generic
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.params.custom_params import global_paths, global_scanning, global_system
from simulation.beam.default_params import beam_params

config = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.sim_tag = "sim_test"
config.scan_tag = "scan_1"

#Options for input_pol_type = ["TQU", "T_only", "noise_only"]
#Default -> "TQU"
config.sim_pol_type = "TQU"

config.add_noise = False

config.do_pencil_beam = True

config.gal_coords = False

config.bolo_config_file = "4_bolos_optimal"

config.pipe_with_map_maker = False

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.t_year = 366*24*60*60.0                    #seconds
config.t_prec = 4*24*60*60.0                     #seconds
config.t_spin = 60.0                             #seconds
config.sampling_rate = 200                                  #Hz

config.alpha = 45.0                                #degrees
config.beta = 45.0                                #degrees

config.oversampling_rate = 1

config.nside_in = 1024

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Noise 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["white", "1_over_f"]
config.noise_type = "white"
config.noise_sigma = 1200                           #\mu K 
config.noise_f_knee = 0
config.noise_alpha = 0

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.bolo_list = ['bolo_0001a', 'bolo_0001b', 'bolo_0002a', 'bolo_0002b']

config.t_segment = 12*60*60.0                      #seconds
config.segment_list = range(8)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Available options for timestream_data_products = ["timestream_data", "scanned_map"]
#"timestream_data" -> will write the simulated pointing, polarisation angle and the simulated signal in the scan data directory
#"scanned_map" -> will write the input hitmap and scanned map in the scan data directory
config.timestream_data_products = ["timestream_data", "scanned_map"]

config.base_dir = global_paths.base_dir 
config.general_data_dir = global_paths.output_dir
config.write_beam = True
config.beam_file_name = "test_beam"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.simulate_beam = True
config.beam_cutoff = 4                                                             #fwhm
config.check_normalisation = True
config.display_beam_settings = True
