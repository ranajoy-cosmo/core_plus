#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import sys

def calculate_params(settings = None):
    if settings is None:
        from settings import settings
    
    if settings.mode is 1:
        settings.t_spin = 360.0*60.0*np.sin(settings.beta_0)*settings.t_sampling/1000.0/settings.theta_co
        settings.t_prec = 360.0*60.0*np.sin(settings.alpha)*settings.t_spin/settings.theta_cross

    if settings.mode is 2:
        settings.theta_cross = 360.0*60.0*np.sin(settings.alpha)*settings.t_spin/settings.t_prec
        settings.theta_co = 360*60*np.sin(settings.beta_0)*settings.t_sampling/1000.0/settings.t_spin

    return settings

def display_params(settings):
    print "alpha : ", np.degrees(settings.alpha), " degrees"
    print "beta : ", np.degrees(settings.beta_0), " degrees"
    print "T flight : ", settings.t_flight/60.0/60.0, "hours"
    print "T precession : ", settings.t_prec/60.0/60.0, "hours"
    print "T spin : ", settings.t_spin, " seconds"
    print "T sampling : ", settings.t_sampling, " milli-seconds"
    print "Scan frequency : ", 1000.0/settings.t_sampling, "Hz"
    print "Theta co : ", settings.theta_co, " arcmin"
    print "Theta cross : ", settings.theta_cross, " arcmin"
    print "Map resolution : ", hp.nside2resol(settings.nside, arcmin = True), " arcmin"

def generate_pointing(settings = None, beta = None):
    if settings is None:
        from settings import settings
    settings = calculate_params(settings)
    if settings.display_params:
        display_params(settings)

    if settings.do_beam_profile_pointing:
        if beta is None:
            print "Beta value is invalid"
            sys.exit()
        else:
            u_init = np.array([np.cos(beta), 0.0, np.sin(beta)])
    else:
        u_init = np.array([np.cos(settings.beta_0), 0.0, np.sin(settings.beta_0)])
    n_steps = int(1000*settings.t_flight/settings.t_sampling)
    print "No. of time steps : ", n_steps
    print "~Memory usage : ", 10*n_steps*8.0/1024/1024/1024, " GB" 
    t_steps = 0.001*settings.t_sampling*np.arange(n_steps)
    w_prec = 2*np.pi/settings.t_prec
    w_spin = 2*np.pi/settings.t_spin
    R = po.Rotation3dOperator("XY'X''", w_prec*t_steps, -1.0*np.full(n_steps, settings.alpha), w_spin*t_steps)
    v = R*u_init
    if settings.write_pointing:
        np.save(settings.output_folder + "pointing_0", v)
    if settings.return_pointing:
        return v

if __name__ == "__main__":
    from default_settings import settings
    settings = calculate_params(settings)
    display_params(settings)
    #generate_pointing(settings)

