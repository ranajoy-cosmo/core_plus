import numpy as np
from simulation.lib.utilities.generic_class import Generic

settings = Generic()

#mode is 1 if we define scan resolution and let scan speed be a function of the resolution
#mode is 2 if we define scan speeds and let the scan resolution be a function of it
settings.mode = 2

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Settings for orientation of spacecraft 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
alpha_deg = 45.0                                    #degrees
settings.alpha = np.deg2rad(alpha_deg)              #radians
beta_deg = 45.0                                     #degrees
settings.beta = np.deg2rad(beta_deg)                #radians

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Settings for time periods of scans
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.t_flight = 10*60*60.0                           #seconds
settings.t_prec = 4*24*60*60.0                     #seconds
settings.t_spin = 60.0                             #seconds
settings.t_sampling = 1.0                          #milli-seconds

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Scan resolutions
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.theta_cross = 10.0                         #arcminutes
settings.theta_co = 10.0                             #arcminutes



