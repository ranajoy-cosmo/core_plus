import numpy as np
import healpy as hp
import os
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.custom_settings import global_paths, global_scanning

settings = Generic()

#mode is 1 if we define scan resolution and let scan speed be a function of the resolution
#mode is 2 if we define scan speeds and let the scan resolution be a function of it
settings.mode = 1 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Settings for orientation of spacecraft 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
alpha_deg = 45.0                                    #degrees
settings.alpha = np.deg2rad(alpha_deg)              #radians
beta_deg = 45.0                                   #degrees
settings.beta = np.deg2rad(beta_deg)            #radians

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Settings for time periods of scans
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.t_year = 365*24*60*60.0                    #seconds
settings.t_flight = 10*60*60.0                       #seconds
settings.t_sampling = 5.0                           #milli-seconds
settings.t_prec = 1*60*60.0                     #seconds
settings.t_spin = 60.0                             #seconds
settings.oversampling_rate = 1

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan resolutions
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.nside = global_scanning.nside_in
settings.theta_cross = hp.nside2resol(settings.nside, arcmin=True)                          #arcmin
settings.theta_co = hp.nside2resol(settings.nside, arcmin=True)                             #arcmin

settings.do_pol = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.output_folder = global_paths.output_folder 
settings.write_pointing = True 
settings.return_pointing = False
settings.display_params = True 
