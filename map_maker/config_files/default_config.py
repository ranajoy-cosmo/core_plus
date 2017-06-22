import healpy as hp
import os
import numpy as np
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.lib.utilities.generic_class import Generic
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.global_config import global_paths

config = Generic()


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.simulate_ts = True/False. True means it will simulate the entire timestream signal and pointing from scratch
#if simulate_ts is True, this is where the timestream data will be written if it is asked for. If False, this is where the timestream data will be read from.
config.simulate_ts = False

config.sim_tag = "sim_test"
config.scan_tag = "scan_test"
config.map_making_tag = "map_recon_test"

#pol_type can be "TQU" OR "QU" OR "T"
config.pol_type = "TQU"

config.take_diff_signal = False
config.subtract_template = False
config.noise_only_map = False

config.notes = "Testing"
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Resolution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.nside_out : nside of the output map
config.nside_out = 1024 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

##Available options for map_making_data_products = ["partial_covariance_matrix", "b_matrix"]
## reconstructed_maps, hitmap, inverse_covariance_maps, covariance_maps are done by default
#config.map_making_data_products = []

config.__dict__.update(global_paths.__dict__)

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

#sim_type can be "signal" OR "gradient"
config.sim_type = "signal"                                                                                                           
#gradient_type can be the following combinations
#"gradient_co" OR "gradient_cross" OR ["gradient_co", "gradient_cross"]
#"gradient_coXgradient_co" OR "gradient_crossXgradient_cross" OR ["gradient_coXgradient_co", "gradient_crossXgradient_cross"]
#"gradient_coXgradient_cross"
#Combinations of the gradients
config.gradient_type = ["gradient_co", "gradient_cross", "gradient_coXgradient_co", "gradient_crossXgradient_cross", "gradient_co_gradi    ent_cross"]

config.add_noise = False

#beam_type can be "pencil" OR "full_simulated" OR "from_file" 
config.beam_type = "pencil"

#coordinate_system can be "ecliptic" OR "galactic"
config.coordinate_system = "ecliptic"

config.bolo_config_file = "simulation.timestream_simulation.bolo_config_files.four_bolos_optimal" 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.t_year = 365.25*24*60*60.0                    #seconds
#* CORE Baseline #*#*#*#*#*
config.t_prec = 4*24*60*60.0                     #seconds
config.t_spin = 120.0                             #seconds
config.sampling_rate = 170                                  #Hz

config.alpha = 30.0                                #degrees
config.beta = 61.0                                #degrees
#*#*#*#*#*#*#*#*#*#*

#* LightBIRD #*#*#*#*#*
#config.t_prec = 93*60.0                     #seconds
#config.t_spin = 600.0                             #seconds
#config.sampling_rate = 10                                  #Hz
#
#config.alpha = 65.0                                #degrees
#config.beta = 30.0                                #degrees
#*#*#*#*#*#*#*#*#*#*

config.oversampling_rate = 1

config.nside_in = 1024

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Noise (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["none", "white_noise", "1_over_f"]
config.noise_type = "none"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Available options for timestream_data_products = ["signal", "pointing_vec", "pol_ang", "noise"]
config.timestream_data_products = ["signal", "pointing_vec", "pol_ang"]
#config.timestream_data_products = config.gradient_type

config.write_beam = True
config.beam_file_name = "test_beam"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters (Only if simulate_ts is True)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#beam_type can be "pencil" OR "full_simulated" OR "from_file" 
config.beam_type = ["pencil"]
config.beam_cutoff = 4                                                             #fwhm
config.check_normalisation = True
