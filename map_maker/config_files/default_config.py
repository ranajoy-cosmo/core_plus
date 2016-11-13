import healpy as hp
import os
import numpy as np
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.lib.utilities.generic_class import Generic
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.params.custom_params import global_paths

config = Generic()


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.simulate_ts = True/False. True means it will simulate the entire timestream signal and pointing from scratch
config.simulate_ts = False
config.map_making_action = "new"

config.sim_tag = "sim_test"
#if simulate_ts is True, this is where the timestream data will be written if it is asked for. If False, this is where the timestream data will be read from.
config.scan_tag = "scan_1"
config.map_making_tag = "map_recon_1"

config.stop_at_inv_cov_map = False

config.take_diff_signal = False
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Resolution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.nside_out : nside of the output map
config.nside_out = 1024 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Available options for map_making_data_products = ["partial_covariance_matrix", "b_matrix"]
# reconstructed_maps, hitmap, inverse_covariance_maps, covariance_maps are done by default
config.map_making_data_products = []

config.base_dir = global_paths.base_dir 
config.general_data_dir = global_paths.output_dir

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#If simulate_ts is true, this is used for simulation as well as map_making, otherwise these are the data segments that are read.

config.bolo_list = ['bolo_0001a', 'bolo_0001b', 'bolo_0002a', 'bolo_0002b']

config.t_segment = 12*60*60.0                      #seconds
config.segment_list = range(8)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for input_pol_type = ["TQU", "T_only", "noise_only"]
#Default -> "TQU"
config.sim_pol_type = "TQU"

config.add_noise = False

config.do_pencil_beam = True

config.gal_coords = False

config.bolo_config_file = "4_bolos_optimal"

config.pipe_with_map_maker = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan (Only if simulate_ts is True)
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
# Noise (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["white_noise", "1_over_f"]
config.noise_type = "white_noise"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Available options for timestream_data_products = ["timestream_data", "scanned_map", "grad_T_scan", "noise"]
#"timestream_data" -> will write the simulated pointing, polarisation angle and the simulated signal in the scan data directory
#"scanned_map" -> will write the input hitmap and scanned map in the scan data directory
config.timestream_data_products = ["timestream_data"]

config.write_beam = True
config.beam_file_name = "test_beam"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.simulate_beam = True
config.beam_cutoff = 4                                                             #fwhm
config.check_normalisation = True
config.display_beam_settings = True
