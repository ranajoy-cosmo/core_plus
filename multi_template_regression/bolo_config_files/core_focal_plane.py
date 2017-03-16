import numpy as np
import os
from simulation.lib.utilities import generic_class

bolo_config = generic_class.Generic()
bolo_config.bolos = {}

scratch_dir = "/scratch1/scratchdirs/banerji"
bolo_config.common_input_map = os.path.join(scratch_dir, "core_maps/r_0001_lensed/sky_map_1024_0.fits")
bolo_config.input_beam_folder = "/global/homes/b/banerji/simulation/beam/nuim"
bolo_config.offset_sigma = 0.0
bolo_config.common_beam_angle = 0.0

bolo_config.common_fwhm = 7.68
bolo_config.common_f_knee = 200                     #mHz
bolo_config.common_white_noise_sigma = 50           #uk*sqrt(s)
bolo_config.common_one_over_f_alpha = 1

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 1a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_1a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 9.875                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 1b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_1b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 9.875                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 1c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_1c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 9.875                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 1d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_1d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 9.875                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 2a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_2a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 19.75                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 2b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_2b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 19.75                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 2c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_2c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 19.75                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 2d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_2d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 19.75                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 3a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_3a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 29.625                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 3b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_3b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 29.625                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 3c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_3c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 29.625                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 3d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_3d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 29.625                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 4a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_4a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 39.5                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 4b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_4b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 39.5                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 4c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_4c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 39.5                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 4d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_4d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 39.5                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 5a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_5a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 49.37                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 5b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_5b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 49.37                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 5c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_5c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 49.37                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 5d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_5d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 49.37                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 6a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_6a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 59.24                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 6b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_6b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 59.24                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 6c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_6c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 59.24                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 6d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_6d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 59.24                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 7a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_7a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 69.12                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 7b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_7b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 69.12                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 7c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_7c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 69.12                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 7d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_7d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 69.12                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 8a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_8a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 78.99                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 8b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_8b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 78.99                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 8c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_8c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 78.99                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 8d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_8d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 78.99                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 9a
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_9a"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 0.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 88.86                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 9b
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_9b"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 90.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 88.86                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 9c
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_9c"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 45.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 88.86                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Bolo 145 9d
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
bolo_name = "bolo_145_9d"
bolo_config.bolos[bolo_name] = generic_class.Generic()

bolo_config.bolos[bolo_name].fwhm_major = bolo_config.common_fwhm                               #arcmins
bolo_config.bolos[bolo_name].fwhm_minor = bolo_config.common_fwhm                                #arcmins
bolo_config.bolos[bolo_name].pol_phase_ini = 135.0                             #degrees
bolo_config.bolos[bolo_name].focal_plane_del_beta = 88.86                        #arcmins 
bolo_config.bolos[bolo_name].beam_angle = bolo_config.common_beam_angle      #degrees

if bolo_config.offset_sigma == 0.0:
    bolo_config.bolos[bolo_name].offset_x = 0.0 
    bolo_config.bolos[bolo_name].offset_y = 0.0
else:
    bolo_config.bolos[bolo_name].offset_x = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)
    bolo_config.bolos[bolo_name].offset_y = np.random.normal(loc=0.0, scale=bolo_config.offset_sigma)

bolo_config.bolos[bolo_name].input_map = bolo_config.common_input_map 

bolo_config.bolos[bolo_name].white_noise_sigma = bolo_config.common_white_noise_sigma
bolo_config.bolos[bolo_name].f_knee = bolo_config.common_f_knee
bolo_config.bolos[bolo_name].one_over_f_alpha = bolo_config.common_one_over_f_alpha
bolo_config.bolos[bolo_name].one_over_f_seed = 1234

bolo_config.bolos[bolo_name].input_beam_file = os.path.join(bolo_config.input_beam_folder, "DRAG_CENTRAL_GRID_linxpol.npy")
