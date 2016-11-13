import numpy as np
import os
from simulation.lib.utilities import generic_class

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

map_folder = "/home/banerji/simulation/maps/bandpass_set_1"
sim_folder = "/home/banerji/simulation/output/orth_test_2"

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

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "group2_map_143.68GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                            #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "group2_map_144.26GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_TEMPLATE"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

#bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "group2_map_350GHz.fits")
bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "sky_diff_TEMPLATE.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001_diff_re
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001_diff_re"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_diff_QU", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE_QU_re
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_TEMPLATE_QU_re"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "rec_TEMPLATE_QU", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
