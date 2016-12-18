import numpy as np
import os
from simulation.lib.utilities import generic_class
from simulation.template_subtraction.config_files.default_config import config

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

map_folder = "/scratch1/scratchdirs/banerji/PSM_OUTPUT/bandpass_mismatch_maps/observations/LiteBird"
scratch_dir = "/scratch1/scratchdirs/banerji/core_output"
sim_folder = os.path.join(scratch_dir, config.sim_tag)

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

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "143.8GHz", "group2_map_143.8GHz.fits")

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

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "144.09GHz", "group2_map_144.09GHz.fits")

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

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "145.63GHz", "group2_map_145.63GHz.fits")

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

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "146.8GHz", "group2_map_146.8GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0003a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0003a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 22.5                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "146.48GHz", "group2_map_146.48GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0003b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0003b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 112.5                            #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "145.17GHz", "group2_map_145.17GHz.fits")

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

#bolo_config.bolos[bolo_name].input_map = os.path.join("/scratch1/scratchdirs/banerji/core_maps/perfect_template_144.fits")
bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "COrE", "350GHz", "group2_map_350GHz.fits")
#bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "COrE_350", "350GHz", "group2_map_350GHz.fits")
#bolo_config.bolos[bolo_name].input_map = os.path.join("/scratch1/scratchdirs/banerji/core_maps/sky_average.fits")

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

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_TEMPLATE_QU", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
