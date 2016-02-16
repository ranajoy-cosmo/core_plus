#! /usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import sys
import os
import importlib
#import h5py
from simulation.lib.quaternion import quaternion

def generate_pointing(scan_params, bolo_params, segment_group, t_start, del_beta=0.0):

    u_view, axis_spin, axis_prec, axis_rev = get_bolo_initial(scan_params, bolo_params, del_beta)

    if del_beta==0.0 and False:
        print "View direction :", u_view
        print "Spin axis :", axis_spin
        print "Precession axis :", axis_prec
        print "Revolution axis :", axis_rev
    print del_beta
    
    n_steps = int(scan_params.t_segment/(scan_params.t_sampling/1000.0))*scan_params.oversampling_rate
    t_steps = t_start + 0.001*(scan_params.t_sampling/scan_params.oversampling_rate)*np.arange(n_steps)

    w_spin = 2*np.pi/scan_params.t_spin
    w_prec = 2*np.pi/scan_params.t_prec
    w_rev = 2*np.pi/scan_params.t_year

    r_spin = quaternion.make_quaternion(w_spin*t_steps, axis_spin)
    r_prec = quaternion.make_quaternion(w_prec*t_steps, axis_prec)
    r_rev = quaternion.make_quaternion(w_rev*t_steps, axis_rev)

    r_total = quaternion.multiply(r_rev, quaternion.multiply(r_prec, r_spin))

    v_view = quaternion.transform(r_total, u_view)

    #if del_beta==0.0:
        #segment_group.create_dataset("vector", data=v_view)
        #np.save(os.path.join(segment_group, "vector"), v_view)

    if scan_params.do_pol:
        pol_init = np.deg2rad(bolo_params.pol_ang)
        pol_ang = ((w_prec + w_spin)*t_steps + pol_init)%np.pi
        #if del_beta==0.0:
            #segment_group.create_dataset("pol_ang", data=pol_ang)
            #np.save(os.path.join(segment_group, "pol_ang"), pol_ang)

    if scan_params.return_pointing:
        if scan_params.do_pol:
            return v_view, pol_ang
        else:
            return v_view


def get_bolo_initial(scan_params, bolo_params, del_beta):
    alpha = np.deg2rad(scan_params.alpha)
    beta = np.deg2rad(scan_params.beta)
    del_x = np.deg2rad(bolo_params.del_x/60.0)
    del_y = np.deg2rad(bolo_params.del_y/60.0)
    del_beta_rad = np.deg2rad(del_beta/60.0)

    u_view = np.array([np.cos(alpha + beta + del_beta_rad), 0.0, np.sin(alpha + beta + del_beta_rad)])
    axis_spin = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
    axis_prec = np.array([1.0, 0.0, 0.0])
    axis_rev = np.array([0.0, 0.0, 1.0])

    return u_view, axis_spin, axis_prec, axis_rev


def calculate_params():
    if scan_params.mode is 1:
        scan_params.t_spin = 360.0*60.0*np.sin(np.radians(scan_params.beta))*scan_params.t_sampling/1000.0/scan_params.theta_co
        scan_params.t_prec = 360.0*60.0*np.sin(np.radians(scan_params.alpha))*scan_params.t_spin/scan_params.theta_cross

    if scan_params.mode is 2:
        scan_params.theta_cross = 360.0*60.0*np.sin(np.radians(scan_params.alpha))*scan_params.t_spin/scan_params.t_prec
        scan_params.theta_co = 360*60*np.sin(np.radians(scan_params.beta))*scan_params.t_sampling/1000.0/scan_params.t_spin

    beam_params.beam_resolution = scan_params.theta_co/scan_params.oversampling_rate


def display_params():
    print "alpha : ", scan_params.alpha, " degrees"
    print "beta : ", scan_params.beta, " degrees"
    print "T flight : ", scan_params.t_flight/60.0/60.0, "hours"
    print "T segment :", scan_params.t_segment/60.0/60.0, "hours"
    print "T precession : ", scan_params.t_prec/60.0/60.0, "hours"
    print "T spin : ", scan_params.t_spin, " seconds"
    print "T sampling : ", scan_params.t_sampling, " milli-seconds"
    print "Scan frequency : ", 1000.0/scan_params.t_sampling, "Hz"
    print "Theta co : ", scan_params.theta_co, " arcmin"
    print "Theta cross : ", scan_params.theta_cross, " arcmin"
    print "Scan resolution for beam integration :", scan_params.theta_co/scan_params.oversampling_rate, "arcmin"
    print "Beam resolution :", beam_params.beam_resolution, "arcmin"
    print "Pixel size for NSIDE =", scan_params.nside, ":", hp.nside2resol(scan_params.nside, arcmin=True), "arcmin"
    n_steps = int(1000*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate
    print "#Samples per segment : ", n_steps
    print "Estimated use of memory : ", 15*n_steps*8.0/1024/1024, "MB"

def run_serial():

    out_dir = os.path.join(scan_params.global_output_dir, "scanning", scan_params.time_stamp)
    #os.makedirs(out_dir)
    
    #root_file = h5py.File(os.path.join(out_dir, "data.hdf5"), libver="latest")

    calculate_params()

    if scan_params.display_params:
        display_params()

    sys.exit()

    for bolo_name in scan_params.bolo_names:
        bolo_group = root_file.create_group(bolo_name)
        bolo_params = importlib.import_module("simulation.timestream_simulation.bolo_params." + bolo_name).bolo_params
        num_segments = int(scan_params.t_flight/scan_params.t_segment)
        print "Bolo Name :", bolo_name
        for segment in range(num_segments):
            segment_name = str(segment+1).zfill(4)
            segment_group = bolo_group.create_group(segment_name)
            t_start = scan_params.t_segment*segment
            print "Segment :", segment_name
            generate_pointing(scan_params, bolo_params, segment_group, t_start)


if __name__ == "__main__":
    from custom_params import scan_params, beam_params
    run_serial()
    #bolo_params = importlib.import_module("simulation.bolo.bolo_params." + "0001").bolo_params
    #scan_params = calculate_params(scan_params, bolo_params)
    #display_params(scan_params, bolo_params)
