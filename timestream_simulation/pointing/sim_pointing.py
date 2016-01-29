#! /usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import sys
import os
import importlib
import h5py
from simulation.lib.quaternion import quaternion

def generate_pointing(settings, bolo_params, segment_group, t_start, del_beta=0.0):
    settings = calculate_params(settings, bolo_params)
    if settings.display_params:
        display_params(settings, bolo_params)

    u_init = get_bolo_initial(bolo_params, del_beta)

    print "Pointing vector in the instrument frame :", int(bolo_params.beta + del_beta/60.0), ((bolo_params.beta + del_beta/60.0)%1)*60

    n_steps = int(settings.t_segment/(settings.t_sampling/1000.0))*settings.oversampling_rate
    t_steps = t_start + 0.001*(settings.t_sampling/settings.oversampling_rate)*np.arange(n_steps)

    w_orbit = 2*np.pi/settings.t_year
    w_prec = 2*np.pi/settings.t_prec
    w_spin = 2*np.pi/settings.t_spin
    alpha = np.deg2rad(bolo_params.alpha)

    r_instrument = quaternion.quaternion_from_euler(w_prec*t_steps, -1.0*alpha, w_spin*t_steps)
    r_orbit = quaternion.make_quaternion(w_orbit*t_steps, np.array([0,0,1]))
    r_total = quaternion.multiply(r_orbit, r_instrument)

    v = quaternion.transform(r_total, u_init)

    if del_beta==0.0:
        segment_group.create_dataset("vector", data=v)

    if settings.do_pol:
        pol_init = np.deg2rad(bolo_params.pol_ang)
        pol_ang = ((w_prec + w_spin)*t_steps + pol_init)%np.pi
        if del_beta==0.0:
            segment_group.create_dataset("pol_ang", data=pol_ang)

    if settings.return_pointing:
        if settings.do_pol:
            return v, pol_ang
        else:
            return v


def get_bolo_initial(bolo_params, del_beta):
    beta = np.deg2rad(bolo_params.beta)
    del_x = np.deg2rad(bolo_params.del_x/60.0)
    del_y = np.deg2rad(bolo_params.del_y/60.0)
    del_beta_rad = np.deg2rad(del_beta/60.0)

    u_init = np.array([np.cos(beta + del_beta_rad), 0.0, np.sin(beta + del_beta_rad)])
    #R = po.Rotation3dOperator('ZY', del_x, del_y)
    #u_init = R*u_init

    return u_init

def calculate_params(settings, bolo_params):
    if settings.mode is 1:
        settings.t_spin = 360.0*60.0*np.sin(bolo_params.beta)*settings.t_sampling/1000.0/settings.theta_co
        settings.t_prec = 360.0*60.0*np.sin(bolo_params.alpha)*settings.t_spin/settings.theta_cross

    if settings.mode is 2:
        settings.theta_cross = 360.0*60.0*np.sin(bolo_params.alpha)*settings.t_spin/settings.t_prec
        settings.theta_co = 360*60*np.sin(bolo_params.beta)*settings.t_sampling/1000.0/settings.t_spin

    return settings

def display_params(settings, bolo_params):
    print "alpha : ", bolo_params.alpha, " degrees"
    print "beta : ", bolo_params.beta, " degrees"
    print "T flight : ", settings.t_flight/60.0/60.0, "hours"
    print "T segment :", settings.t_segment/60.0/60.0, "hours"
    print "T precession : ", settings.t_prec/60.0/60.0, "hours"
    print "T spin : ", settings.t_spin, " seconds"
    print "T sampling : ", settings.t_sampling, " milli-seconds"
    print "Scan frequency : ", 1000.0/settings.t_sampling, "Hz"
    print "Theta co : ", settings.theta_co, " arcmin"
    print "Theta cross : ", settings.theta_cross, " arcmin"
    n_steps = int(1000*settings.t_segment/settings.t_sampling)*settings.oversampling_rate
    print "#Samples per segment : ", n_steps
    print "Estimated use of memory : ", 15*n_steps*8.0/1024/1024, "MB"

def run_serial(settings):

    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp)
    os.makedirs(out_dir)
    
    root_file = h5py.File(os.path.join(out_dir, "data.hdf5"), libver="latest")

    for bolo_name in settings.bolo_names:
        bolo_group = root_file.create_group(bolo_name)
        bolo_params = importlib.import_module("simulation.timestream_simulation.bolo.bolo_params." + bolo_name).bolo_params
        num_segments = int(settings.t_flight/settings.t_segment)
        print "Bolo Name :", bolo_name

        for segment in range(num_segments):
            segment_name = str(segment+1).zfill(4)
            segment_group = bolo_group.create_group(segment_name)
            t_start = settings.t_segment*segment
            print "Segment :", segment_name
            generate_pointing(settings, bolo_params, segment_group, t_start)


if __name__ == "__main__":
    from custom_settings import settings
    run_serial(settings)
    #bolo_params = importlib.import_module("simulation.bolo.bolo_params." + "0001").bolo_params
    #settings = calculate_params(settings, bolo_params)
    #display_params(settings, bolo_params)
