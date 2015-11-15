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
    print "Scan speed : ", 360*np.sin(settings.beta)/settings.t_spin, " degrees/sec"
    print "Scan frequency : ", 1000.0/settings.t_sampling, "Hz"
    print "Theta co : ", settings.theta_co, " arcmin"
    print "Theta cross : ", settings.theta_cross, " arcmin"

def generate_pointing(settings=None, del_beta=0):
    if settings is None:
        from custom_settings import settings
    settings = calculate_params(settings)
    if settings.display_params:
        display_params(settings)

    u_init = np.array([np.cos(settings.beta + del_beta), 0.0, np.sin(settings.beta + del_beta)])
    print int(np.degrees(settings.beta + del_beta)), (np.degrees(settings.beta + del_beta)%1)*60

    n_steps = int(1000*settings.t_flight/settings.t_sampling)*settings.oversampling_rate
    t_steps = 0.001*settings.t_sampling*np.arange(n_steps)/settings.oversampling_rate
    w_orbit = 2*np.pi/settings.t_year
    w_prec = 2*np.pi/settings.t_prec
    w_spin = 2*np.pi/settings.t_spin
    R = po.Rotation3dOperator("XY'X''", -1.0*w_prec*t_steps, -1.0*np.full(n_steps, settings.alpha), w_spin*t_steps)
    v = R*u_init
    R = po.Rotation3dOperator("Z", w_orbit*t_steps)
    v = R*v
    #lat = np.pi/2 + np.random.random(n_steps)*np.pi*10/180
    #lon = np.random.random(n_steps)*np.pi*10/180
    #v = hp.ang2vec(lat, lon)

    if settings.do_pol is True and del_beta == 0.0:
        #pol_ang = (w_prec + w_spin)*t_steps%np.pi
        pol_ang = np.random.random(n_steps)*np.pi
        if settings.write_pointing:
            np.save(os.path.join(settings.output_folder, "pointing.npy"), v[::settings.oversampling_rate])
            np.save(os.path.join(settings.output_folder, "pol_angle.npy"), pol_ang[::settings.oversampling_rate])
            np.save(os.path.join(settings.output_folder, "times.npy"), t_steps[::settings.oversampling_rate])
        if settings.return_pointing:
            return v, pol_ang

    else:
        if settings.write_pointing and del_beta == 0.0:
            np.save(os.path.join(settings.output_folder, "pointing.npy"), v[::settings.oversampling_rate])
            np.save(os.path.join(settings.output_folder, "times.npy"), t_steps[::settings.oversampling_rate])
        if settings.return_pointing:
            return v

if __name__ == "__main__":
    from custom_settings import settings
    generate_pointing(settings)

