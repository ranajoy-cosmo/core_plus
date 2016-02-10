#!/usr/bin/env python 

from default_params import *

scan_params.bolo_names = ['bolo_0001', 'bolo_0002', 'bolo_0003', 'bolo_0004']

scan_params.t_flight = 16*24*60*60.0
scan_params.t_segment = 8*60*60.0

scan_params.input_map = os.path.join(global_paths.maps_dir, "sky_map_4096_0_pixwin_T.fits") 

beam_params.do_pencil_beam = False

scan_params.oversampling_rate = 1 
