import numpy as np
import os
from simulation.lib.utilities import generic_class
from simulation.multi_template_regression.config_files.default_config import config

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

scratch_dir = "/scratch1/scratchdirs/banerji"
map_folder = os.path.join(scratch_dir, "PSM_OUTPUT/bandpass_mismatch_maps/observations/CORE_new")
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
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "144.17GHz", "group3_map_144.17GHz.fits")

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
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "145.85GHz", "group3_map_145.85GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE TD
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001_TEMPLATE_TD"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "340GHz", "group1_map_340GHz.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE_SYNC
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001_TEMPLATE_SYNC"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "90GHz", "group2_map_90GHz.fits")

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
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_diff_QU_1", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE TD QU rescan
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001_TEMPLATE_TD_re"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_TEMPLATE_TD_1_QU", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo TEMPLATE SYNC QU rescan
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001_TEMPLATE_SYNC_re"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 0.0                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 0.0                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees

bolo_config.bolos[bolo_name].offset_x = 0.0 
bolo_config.bolos[bolo_name].offset_y = 0.0
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0

bolo_config.bolos[bolo_name].input_map = os.path.join(sim_folder, "rec_TEMPLATE_SYNC_1_QU", "sky_map.fits")

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
