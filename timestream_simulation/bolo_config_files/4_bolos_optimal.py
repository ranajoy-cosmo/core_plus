import numpy as np
import os
from simulation.lib.utilities import generic_class

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

bolo_config.common_input_map = "/home/banerji/simulation/maps/r_0001/sky_map_1024_0.fits"
map_folder = "/home/banerji/simulation/maps/bandpass_set_1"
bolo_config.input_beam_folder = "/project/projectdirs/core/sim/beams/planck_7.68_arcmin/square"
bolo_config.common_input_beam = os.path.join(bolo_config.input_beam_folder, "beam_217-5a_uv_rescaled_fwhm_7.68_arcmin.fits")
bolo_config.offset_sigma = 0.0
bolo_config.common_beam_angle = 0.0

bolo_config.common_f_knee = 200                     #mHz
bolo_config.common_white_noise_sigma = 50           #uk*sqrt(s)
bolo_config.common_one_over_f_alpha = 1

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 
#bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "group2_map_143.42GHz.fits")
bolo_config.bolos[bolo_name].input_map = "/home/banerji/simulation/output/bandpass_test_gal/rec_diff_QU/sky_map.fits" 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = 0.0 
#bolo_config.bolos[bolo_name].offset_y = 0.0

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "beam_217-5a_uv_rescaled_fwhm_7.68_arcmin.fits")

bolo_config.bolos[bolo_name].beam_angle = 0.0                                  #degrees

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                            #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "group2_map_146.88GHz.fits")
#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 2345

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
#bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "beam_217-5b_uv_rescaled_fwhm_7.68_arcmin.fits")

bolo_config.bolos[bolo_name].beam_angle = 0.0      #degrees

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0002a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                            #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 
#bolo_config.bolos[bolo_name].input_beam_file = bolo_config.common_input_beam

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 3456

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
#bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "beam_217-6a_uv_rescaled_fwhm_7.68_arcmin.fits")

bolo_config.bolos[bolo_name].beam_angle = 0.0                         #degrees

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0002b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                           #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
#bolo_config.bolos[bolo_name].input_beam_file = bolo_config.common_input_beam

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 4567

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
#bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "beam_217-6b_uv_rescaled_fwhm_7.68_arcmin.fits")

bolo_config.bolos[bolo_name].beam_angle = 0.0      #degrees

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 350
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_350"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 7.7266                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                           #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = os.path.join(map_folder, "group2_map_350GHz_T_only.fits")
bolo_config.bolos[bolo_name].input_map = "/home/banerji/simulation/output/bandpass_test_gal/rec_350_QU/sky_map.fits" 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 4567

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
#bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "beam_217-6b_uv_rescaled_fwhm_7.68_arcmin.fits")

bolo_config.bolos[bolo_name].beam_angle = 0.0      #degrees

