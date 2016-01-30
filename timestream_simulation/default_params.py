import healpy as hp
import os
from simulation.lib.utilities.generic_class import Generic
from simulation.lib.utilities.time_util import get_time_stamp
from simulation.settings.custom_settings import global_paths, global_scanning, global_system
from simulation.beam.default_params import beam_params

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

scan_params.bolo_names = ['bolo_0001']
scan_params.nside = global_scanning.nside_in
scan_params.do_only_T = False
scan_params.oversampling_rate = 1

scan_params.mode = 2 
scan_params.t_year = 365*24*60*60.0                    #seconds
scan_params.t_flight = 4*24*60*60.0                       #seconds
scan_params.t_segment = 12*60*60.0
scan_params.t_sampling = 5.0                           #milli-seconds
scan_params.t_prec = 4*24*60*60.0                     #seconds
scan_params.t_spin = 60.0                             #seconds

scan_params.theta_cross = hp.nside2resol(scan_params.nside, arcmin=True)                          #arcmin
scan_params.theta_co = hp.nside2resol(scan_params.nside, arcmin=True)                             #arcmin

scan_params.do_pol = True
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

scan_params.time_stamp = get_time_stamp()
scan_params.base_dir = global_paths.base_dir 
scan_params.global_output_dir = global_paths.output_dir
scan_params.input_map = os.path.join(global_paths.maps_dir, "sky_map_4096_0_pixwin.fits")
scan_params.bolo_param_dir = os.path.join(global_paths.base_dir, "bolo", "bolo_params")
scan_params.time_stamp = get_time_stamp()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

beam_params.do_pencil_beam = True
