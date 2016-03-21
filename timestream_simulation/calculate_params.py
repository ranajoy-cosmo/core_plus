#!/usr/bin/env python 

import numpy as np
import healpy as hp
from simulation.lib.geometry.conversions import *
from custom_params import *

scan_params.theta_co = scan_params.beam_fwhm/scan_params.samples_per_beam                   #arcmin
scan_radius = 360.0*60*np.sin(np.radians(scan_params.beta))                                 #arcmin
scan_speed = scan_radius/scan_params.t_spin                                                 #arcmin/seconds
scan_params.sampling_rate = scan_speed/scan_params.theta_co

if scan_params.mode is 1:
    scan_velocity = scan_params.theta_co*scan_params.sampling_rate                                             #arcmin/second
    scan_params.t_spin = 360.0*60.0*np.sin(np.radians(scan_params.beta))*scan_params.t_sampling/1000.0/scan_params.theta_co
    scan_params.t_prec = 360.0*60.0*np.sin(np.radians(scan_params.alpha))*scan_params.t_spin/scan_params.theta_cross

if scan_params.mode is 2:
    scan_velocity = 360.0*60.0*sin(np.radians(scan_params.beta))/scan_params.t_spin                                     #arcmin/second
    scan_params.theta_co = scan_velocity/scan_params.sampling_rate                                                      #arcmin
    scan_params.theta_cross = 360.0*60.0*np.sin(np.radians(scan_params.alpha))*scan_params.t_spin/scan_params.t_prec    #arcmin

beam_params.beam_resolution = scan_params.theta_co/scan_params.oversampling_rate
scan_params.scan_resolution = scan_params.theta_co/scan_params.oversampling_rate
pix_width = np.degrees(hp.nside2resol(scan_params.nside))
scan_speed = 360.0*np.sin(np.radians(scan_params.beta))/scan_params.t_spin
scan_params.filter_cutoff = scan_speed/pix_width

def display_params():
    print "alpha : ", scan_params.alpha, " degrees"
    print "beta : ", scan_params.beta, " degrees"
    print "Mode :", scan_params.mode
    print "T flight : ", scan_params.t_flight/60.0/60.0, "hours /", scan_params.t_flight/60.0/60.0/24, "days"
    print "T segment :", scan_params.t_segment/60.0/60.0, "hours /", scan_params.t_segment/60.0/60.0/24, "days"
    print "T precession : ", scan_params.t_prec/60.0/60.0, "hours"
    print "T spin : ", scan_params.t_spin, " seconds"
    print "T sampling : ", scan_params.t_sampling, " milli-seconds"
    print "Scan frequency : ", 1000.0/scan_params.t_sampling, "Hz"
    print "Theta co : ", scan_params.theta_co, " arcmin"
    print "Theta cross : ", scan_params.theta_cross, " arcmin"
    print "Scan resolution for beam integration :", scan_params.scan_resolution, "arcmin"
    print "Beam resolution :", beam_params.beam_resolution, "arcmin"
    print "Pixel size for NSIDE =", scan_params.nside, ":", hp.nside2resol(scan_params.nside, arcmin=True), "arcmin"
    print "Filter cutoff :", scan_params.filter_cutoff
    n_steps = int(1000*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate
    print "#Samples per segment : ", n_steps
    print "Size of signal array : ", n_steps*8.0/1024/1024, "MB"
    print "Estimated use of memory : ", 15*n_steps*8.0/1024/1024, "MB"


display_params()
