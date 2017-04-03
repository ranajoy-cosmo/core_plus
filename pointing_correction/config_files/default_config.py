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

config.sim_tag = "pointing_correction_test_4"
config.scan_tag = "scan"
config.resim = False

#coordinate_system can be "ecliptic" OR "galactic"
config.coordinate_system = "ecliptic"

config.bolo_config_file = "simulation.pointing_correction.bolo_config_files.bolo_config" 

config.overwrite = False

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.t_year = 365.25*24*60*60.0                    #seconds

#* CORE Baseline #*#*#*#*#*
config.scan_strategy_name = "CORE Baseline"
config.t_prec = 4*24*60*60.0                     #seconds
config.t_spin = 120.0                             #seconds
config.sampling_rate = 170                                  #Hz

config.alpha = 30.0                                #degrees
config.beta = 65.0                                #degrees
#*#*#*#*#*#*#*#*#*#*

config.oversampling_rate = 1

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Resolution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
config.nside_in = 4096
config.nside_out = 4096

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Noise 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["none", "white", "1_over_f"]
config.noise_type = "none"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#config.t_segment = 1.90234375*24*60*60.0                      #seconds
#config.segment_list = range(192)
config.t_segment = 0.4755859375*24*60*60.0                      #seconds
config.segment_list = range(768)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.__dict__.update(global_paths.__dict__)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Input Maps 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

input_map = os.path.join(config.maps_dir, "r_0001_lensed/sky_map_4096_0.fits")
diff_rescan_map = os.path.join(config.general_output_dir, config.sim_tag, "rec_diff", "sky_map.fits")
tm_grad_co_rescan_map = os.path.join(config.general_output_dir, config.sim_tag, "tm_grad_co", "sky_map.fits")
tm_grad_cross_rescan_map = os.path.join(config.general_output_dir, config.sim_tag, "tm_grad_cross", "sky_map.fits")

config.fill_empty_pixels = True

config.beam_cutoff = 4                                                             #fwhm
