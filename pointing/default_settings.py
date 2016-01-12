import numpy as np
import healpy as hp
import os
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.custom_settings import global_paths, global_scanning
from simulation.lib.utilities.time_util import get_time_stamp

settings = Generic()

#mode is 1 if we define scan resolution and let scan speed be a function of the resolution
#mode is 2 if we define scan speeds and let the scan resolution be a function of it
settings.mode = 2 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Settings for time periods of scans
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.t_year = 365*24*60*60.0                    #seconds
settings.t_flight = 4*24*60*60.0                       #seconds
settings.t_segment = 12*60*60.0
settings.t_sampling = 5.0                           #milli-seconds
settings.t_prec = 4*24*60*60.0                     #seconds
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
settings.bolo_names = ['0001', '0002']
settings.time_stamp = get_time_stamp()
settings.global_output_dir = global_paths.output_dir
settings.write_pointing = True 
settings.return_pointing = False
settings.display_params = False 
