import numpy as np
import os
from simulation.lib.utilities import generic_class

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

map_folder = "/home/banerji/simulation/output/bandpass_test_gal_6"

bolo_config.common_white_noise_sigma = 50           #uk*sqrt(s)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "rec_diff", "sky_map.fits") 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 350
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_TEMPLATE"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "rec_TEMPLATE_QU", "sky_map.fits") 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
