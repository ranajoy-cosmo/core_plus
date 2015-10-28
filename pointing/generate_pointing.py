#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import sys, os

def calculate_params(settings):
    if settings.mode is 1:
        settings.t_spin = 360.0*60.0*np.sin(settings.beta)*settings.t_sampling/1000.0/settings.theta_co
        settings.t_prec = 360.0*60.0*np.sin(settings.alpha)*settings.t_spin/settings.theta_cross

    if settings.mode is 2:
        settings.theta_cross = 360.0*60.0*np.sin(settings.alpha)*settings.t_spin/settings.t_prec
        settings.theta_co = 360*60*np.sin(settings.beta)*settings.t_sampling/1000.0/settings.t_spin

    return settings

def display_params(settings):
    print "alpha : ", np.degrees(settings.alpha), " degrees"
    print "beta : ", np.degrees(settings.beta), " degrees"
    print "T flight : ", settings.t_flight/60.0/60.0, "hours"
    print "T precession : ", settings.t_prec/60.0/60.0, "hours"
    print "T spin : ", settings.t_spin, " seconds"
    print "T sampling : ", settings.t_sampling, " milli-seconds"
    print "Scan frequency : ", 1000.0/settings.t_sampling, "Hz"
    print "Theta co : ", settings.theta_co, " arcmin"
    print "Theta cross : ", settings.theta_cross, " arcmin"

def generate_pointing(settings=None, del_beta=0):
    if settings is None:
        from local_settings import settings
    settings = calculate_params(settings)
    if settings.display_params:
        display_params(settings)

    u_init = np.array([np.cos(settings.beta + del_beta), 0.0, np.sin(settings.beta + del_beta)])
    print np.degrees(settings.beta +  del_beta)

    n_steps = int(1000*settings.t_flight/settings.t_sampling)
    t_steps = 0.001*settings.t_sampling*np.arange(n_steps)
    w_prec = 2*np.pi/settings.t_prec
    w_spin = 2*np.pi/settings.t_spin
    R = po.Rotation3dOperator("XY'X''", -1.0*w_prec*t_steps, -1.0*np.full(n_steps, settings.alpha), w_spin*t_steps)
    v = R*u_init

    if settings.do_pol is True and del_beta == 0.0:
        pol_ang = (w_prec + w_spin)*t_steps%np.pi
        print "Flag 1"
        if settings.write_pointing:
            print "Flag 2"
            np.save(os.path.join(settings.output_folder, "pointing.npy"), v)
            np.save(os.path.join(settings.output_folder, "pol_angle.npy"), pol_ang)
            np.save(os.path.join(settings.output_folder, "times.npy"), t_steps)
        if settings.return_pointing:
            print "Flag 3"
            return v, pol_ang

    else:
        print "Flag 4"
        if settings.write_pointing and del_beta == 0.0:
            print "Flag 5"
            np.save(os.path.join(settings.output_folder, "pointing.npy"), v)
            np.save(os.path.join(settings.output_folder, "times.npy"), t_steps)
        if settings.return_pointing:
            print "Flag 6"
            return v

if __name__ == "__main__":
    from local_settings import settings
    generate_pointing(settings)
