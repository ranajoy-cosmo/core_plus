import numpy as np
import os
from simulation.lib.utilities import generic_class
from simulation.template_subtraction.config_files.default_config import config

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

scratch_dir = "/scratch1/scratchdirs/banerji"
map_folder = os.path.join(scratch_dir, "PSM_OUTPUT/bandpass_mismatch_maps/observations/CORE")
sim_folder = os.path.join(scratch_dir, "core_output", config.sim_tag)

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

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "143.677GHz", "group1_map_143.677GHz.fits")

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

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "146.041GHz", "group1_map_146.041GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0002a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "146.061GHz", "group1_map_146.061GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0002b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                            #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "146.334GHz", "group1_map_146.334GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE 1
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001_TEMPLATE"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "350GHz", "group2_map_350GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE 2
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002_TEMPLATE"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "350GHz", "group2_map_350GHz.fits")

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

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_diff_QU_1", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0002_diff_re
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002_diff_re"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_diff_QU_2", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE_QU_re_1
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001_TEMPLATE_QU_re"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_TEMPLATE_QU_1", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE_QU_re_2
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002_TEMPLATE_QU_re"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_TEMPLATE_QU_2", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
