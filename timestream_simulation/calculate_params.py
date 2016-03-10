#!/usr/bin/env python 

import numpy as np
import healpy as hp
from custom_params import *

if scan_params.mode is 1:
    scan_params.t_spin = 360.0*60.0*np.sin(np.radians(scan_params.beta))*scan_params.t_sampling/1000.0/scan_params.theta_co
    scan_params.t_prec = 360.0*60.0*np.sin(np.radians(scan_params.alpha))*scan_params.t_spin/scan_params.theta_cross

if scan_params.mode is 2:
    scan_params.theta_cross = 360.0*60.0*np.sin(np.radians(scan_params.alpha))*scan_params.t_spin/scan_params.t_prec
    scan_params.theta_co = 360*60*np.sin(np.radians(scan_params.beta))*scan_params.t_sampling/1000.0/scan_params.t_spin

beam_params.beam_resolution = scan_params.theta_co/scan_params.oversampling_rate
scan_params.scan_resolution = scan_params.theta_co/scan_params.oversampling_rate
pix_width = np.degrees(hp.nside2resol(scan_params.nside))
scan_speed = 360.0*np.sin(np.radians(scan_params.beta))/scan_params.t_spin
scan_params.filter_cutoff = scan_speed/pix_width

