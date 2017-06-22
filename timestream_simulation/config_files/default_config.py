import healpy as hp
import os
import numpy as np
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.lib.utilities.generic_class import Generic
from simulation.global_config import global_paths

config = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.simulate_ts = True
config.sim_tag = "sim_test"
config.scan_tag = "scan_test"

#Options for input_pol_type = ["TQU", "QU", "T", "_QU", "noise_only"]
#Default -> "TQU"
# A _ means that the input map has a component which will not be read"
config.sim_pol_type = "TQU"

#coordinate_system can be "ecliptic" OR "galactic"
config.coordinate_system = "ecliptic"

config.bolo_config_file = "simulation.timestream_simulation.bolo_config_files.four_bolos_optimal"

#sim_type can be "signal" OR "gradient"
config.sim_type = ["signal"]
#if sim_type is chosen as "gradient", gradient_type can be the following combinations
#"gradient_co" OR "gradient_cross" OR ["gradient_co", "gradient_cross"]
#"gradient_coXgradient_co" OR "gradient_crossXgradient_cross" OR ["gradient_coXgradient_co", "gradient_crossXgradient_cross"]
#"gradient_coXgradient_cross"
#Combinations of the gradients
config.gradient_type = []

config.overwrite = False

config.notes = "Testing"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.t_year = 365.25*24*60*60.0                    #seconds
#* CORE Baseline #*#*#*#*#*
#config.t_prec = 4*24*60*60.0                     #seconds
#config.t_spin = 120.0                             #seconds
#config.sampling_rate = 170                                  #Hz
#
#config.alpha = 30.0                                #degrees
#config.beta = 61.0                                #degrees
#*#*#*#*#*#*#*#*#*#*

#* LightBIRD #*#*#*#*#*
config.t_prec = 93*60.0                     #seconds
config.t_spin = 600.0                             #seconds
config.sampling_rate = 10                                  #Hz

config.alpha = 65.0                                #degrees
config.beta = 30.0                                #degrees
#*#*#*#*#*#*#*#*#*#*

config.oversampling_rate = 1

config.nside_in = 1024

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Noise 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["none", "white", "1_over_f"]
config.noise_type = "white"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.bolo_list = ['bolo_0001a', 'bolo_0001b', 'bolo_0002a', 'bolo_0002b']

config.t_segment = 12*60*60.0                      #seconds
config.segment_list = range(8)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Available options for timestream_data_products = ["signal", "pointing_vec", "pol_ang", "noise"]
config.timestream_data_products = ["signal", "pointing_vec", "pol_ang"]

config.__dict__.update(global_paths.__dict__)

config.write_beam = True
config.beam_file_name = "test_beam"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#beam_type can be "pencil" OR "full_simulated" OR "from_file" 
config.beam_type = ["pencil"]
config.beam_cutoff = 4                                                             #fwhm
config.check_normalisation = True
