#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import sys, os, importlib
#import simulation.lib.utilities.memory_management as mem

machine='edison'
num_proc = 24

def generate_pointing(segment, settings, bolo_params, del_beta=0.0):
    #mem.check_memory("Pointing flag 1", macnine, num_proc) 
    settings = calculate_params(settings, bolo_params)
    if settings.display_params:
        display_params(settings, bolo_params)

    u_init = get_bolo_initial(bolo_params, del_beta)

    print int(bolo_params.beta + del_beta/60.0), ((bolo_params.beta + del_beta/60.0)%1)*60

    n_steps = int(1000*settings.t_segment/settings.t_sampling)*settings.oversampling_rate
    t_start = settings.t_segment*segment
    t_steps = t_start + 0.001*(settings.t_sampling/settings.oversampling_rate)*np.arange(n_steps)

    print "Time range : [", t_steps[0], ", ", t_steps[-1], "]" 

    w_orbit = 2*np.pi/settings.t_year
    w_prec = 2*np.pi/settings.t_prec
    w_spin = 2*np.pi/settings.t_spin
    alpha = np.deg2rad(bolo_params.alpha)
    #mem.check_memory("Pointing flag 2", macnine, num_proc) 

    R = po.Rotation3dOperator("XY'X''", -1.0*w_prec*t_steps, -1.0*np.full(n_steps, alpha), w_spin*t_steps)
    v = R*u_init
    #mem.check_memory("Pointing flag 3", macnine, num_proc) 
    R = po.Rotation3dOperator("Z", w_orbit*t_steps)
    v = R*v
    #mem.check_memory("Pointing flag 4", macnine, num_proc) 
    #lat = np.pi/2 + np.random.random(n_steps)*np.pi*10/180
    #lon = np.random.random(n_steps)*np.pi*10/180
    #v = hp.ang2vec(lat, lon)

    if settings.do_pol:
        pol_init = np.deg2rad(bolo_params.pol_ang)
        pol_ang = ((w_prec + w_spin)*t_steps + pol_init)%np.pi  
        #pol_ang = np.random.random(n_steps)*np.pi

    if settings.write_pointing and del_beta==0.0:
        if settings.do_pol: 
            write_pointing(settings, bolo_params, segment, v, pol_ang)
        else:
            write_pointing(settings, bolo_params, segment, v)

    #mem.check_memory("Pointing flag 5", macnine, num_proc) 
    if settings.return_pointing:
        if settings.do_pol:
            return v, pol_ang
        else:
            return v

def make_output_dirs(settings):
    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp)
    for bolo_name in settings.bolo_names:
        os.makedirs(os.path.join(out_dir, bolo_name))

def write_pointing(settings, bolo_params, segment, v, pol=None):
    out_dir = os.path.join(settings.global_output_dir, "scanning", settings.time_stamp, bolo_params.bolo_name) 
    out_file = str(segment+1).zfill(4)
    v_file = os.path.join(out_dir, 'vec_' + out_file)
    np.save(v_file, v[::settings.oversampling_rate])

    if settings.do_pol:
        pol_file = os.path.join(out_dir, 'pol_' + out_file)
        np.save(pol_file, pol[::settings.oversampling_rate])

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
    make_output_dirs(settings)
    for bolo_name in settings.bolo_names:
        bolo_params = importlib.import_module("simulation.bolo.bolo_params." + bolo_name).bolo_params 
        num_segments = int(settings.t_flight/settings.t_segment)
        print "Doing bolo : ", bolo_name
        for segment in range(num_segments):
            print "Segment : ", segment
            generate_pointing(segment, settings, bolo_params) 

if __name__ == "__main__":
    from custom_settings import settings
    #run_serial(settings)
    bolo_params = importlib.import_module("simulation.bolo.bolo_params." + "0001").bolo_params 
    settings = calculate_params(settings, bolo_params)
    display_params(settings, bolo_params)

