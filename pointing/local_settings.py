import numpy as np
import healpy as hp
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.global_settings import global_settings

settings = Generic()

#mode is 1 if we define scan resolution and let scan speed be a function of the resolution
#mode is 2 if we define scan speeds and let the scan resolution be a function of it
settings.mode = 2 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Settings for orientation of spacecraft 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
alpha_deg = 35.0                                    #degrees
settings.alpha = np.deg2rad(alpha_deg)              #radians
beta_deg_0 = 55.0                                   #degrees
settings.beta_0 = np.deg2rad(beta_deg_0)            #radians

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Settings for time periods of scans
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.t_flight = 1*60*60.0                       #seconds
settings.t_sampling = 5.0                           #milli-seconds
settings.t_prec = 4*24*60*60.0                     #seconds
settings.t_spin = 60.0                             #seconds

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan resolutions
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.nside = 2048
settings.theta_cross = 1.0                          #arcmin
settings.theta_co = 0.6                             #arcmin

settings.do_beam_profile_pointing = False
settings.do_polarisation_modulation = False

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
tag = "test"
base_folder = global_settings.base_folder
settings.output_folder = global_settings.output_folder + tag + '/'
settings.write_pointing = True
settings.return_pointing = False
settings.display_params = True