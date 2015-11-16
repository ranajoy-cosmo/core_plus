import numpy as np
import healpy as hp
import os
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.custom_settings import global_paths, global_scanning
import simulation.pointing.default_settings as pointing_settings
import simulation.beam.default_settings as beam_settings
from simulation.lib.utilities.time_util import get_time_stamp

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Parameters 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

settings.bolo_names = ['0001']
settings.nside_in = global_scanning.nside_in
settings.load_pointing = True
settings.do_pol = True
settings.oversampling_rate = 1

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

settings.base_dir = global_paths.base_dir 
settings.global_output_dir = global_paths.output_dir
settings.input_map = os.path.join(global_paths.maps_dir, "sky_map_4096_0.fits")
settings.bolo_param_dir = os.path.join(global_paths.base_dir, "bolo", "bolo_params")
settings.time_stamp = get_time_stamp() 

settings.write_signal = True 
settings.write_scanned_map = True
settings.write_hitmap = True
settings.pipe_with_map_maker = False 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Pointing Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

pointing_params = pointing_settings.settings
pointing_params.t_flight = 200*60*60.0                       #seconds
pointing_params.write_pointing = True
pointing_params.return_pointing = True
pointing_params.display_params = False
pointing_params.do_pol = True
pointing_params.oversampling_rate = settings.oversampling_rate
pointing_params.time_stamp = settings.time_stamp

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

beam_params = beam_settings.settings 
beam_params.do_pencil_beam = True
beam_params.beam_resolution = pointing_params.theta_co/settings.oversampling_rate 
