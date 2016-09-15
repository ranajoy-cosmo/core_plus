import numpy as np
from simulation.lib.utilities import generic_class

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

bolo_config.common_input_map = "/global/homes/b/banerji/simulation/maps/r_001_lensed/sky_map_4096_0.fits"
#bolo_config.common_input_map = "/Users/banerji/CORE+/simulation/maps/r_001/sky_map_1024_8.fits"
bolo_config.common_input_beam = "/global/homes/b/banerji/grasp_beams/square/beam_217-5a_uv_rescaled_fwhm_5.79_arcmin.npy"
bolo_config.offset_sigma = 0.0
bolo_config.common_beam_angle = 0.0

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 5.856                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 5.856                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = bolo_config.common_input_beam

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = 0.0 
#bolo_config.bolos[bolo_name].offset_y = 0.0

#bolo_config.bolos[bolo_name].input_map = 
bolo_config.bolos[bolo_name].input_beam_file = "/global/homes/b/banerji/grasp_beams/square/beam_217-5a_uv_rescaled_fwhm_5.79_arcmin.npy" 

#bolo_config.bolos[bolo_name].beam_angle =                                   #degrees

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0001b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0001b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 5.856                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 5.856                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                            #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = bolo_config.common_input_beam

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
#bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = "/global/homes/b/banerji/grasp_beams/square/beam_217-5b_uv_rescaled_fwhm_5.79_arcmin.npy" 

#bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0002a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 5.856                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 5.856                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                            #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = bolo_config.common_input_beam

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
#bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = "/global/homes/b/banerji/grasp_beams/square/beam_217-6a_uv_rescaled_fwhm_5.79_arcmin.npy" 

#bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 0002b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_0002b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = 5.856                                #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = 5.856                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                           #degrees
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = bolo_config.common_input_beam

#Uncomment and change the next if this bolo needs specific values different from the common values
#bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
#bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

#bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map
bolo_config.bolos[bolo_name].input_beam_file = "/global/homes/b/banerji/grasp_beams/square/beam_217-6b_uv_rescaled_fwhm_5.79_arcmin.npy" 

#bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

