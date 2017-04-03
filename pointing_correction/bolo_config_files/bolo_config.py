import numpy as np
import os
from simulation.lib.utilities import generic_class
from simulation.pointing_correction.config_files.default_config import config

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

scratch_dir = "/scratch1/scratchdirs/banerji"
sim_dir = os.path.join(scratch_dir, "core_output", config.sim_tag)

bolo_config.offset_sigma = 0.0
bolo_config.common_beam_angle = 0.0

bolo_config.common_f_knee = 200                     #mHz
bolo_config.common_white_noise_sigma = 39.9/2.0           #uk*sqrt(s)
bolo_config.common_one_over_f_alpha = 1

bolo_config.common_focal_plane_del_beta = 0.0       #arcmins
bolo_config.common_fwhm = 8.0                       #arcmins

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].name = bolo_name                                                   #arcmins
bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_effective = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = 0.0      #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0 

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

#Uncomment and change the next if this bolo needs specific values different from the common values
bolo_config.bolos[bolo_name].offset_x = -7.8
bolo_config.bolos[bolo_name].offset_y = 13.2

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].name = bolo_name                                                   #arcmins
bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_effective = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                            #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = 0.0                                 #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 0.0 

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 2345

#Uncomment and change the next if this bolo needs specific values different from the common values
bolo_config.bolos[bolo_name].offset_x = 5.1
bolo_config.bolos[bolo_name].offset_y = 2.5
