#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
import copy
import os
import shutil
import importlib
import time
from memory_profiler import profile
from simulation.lib.quaternion import quaternion
from pysimulators import ProjectionOperator, BeamGaussian
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
from simulation.beam.beam_kernel import get_beam
from simulation.lib.utilities.time_util import get_time_stamp
import simulation.lib.numericals.filters as filters
import simulation.lib.data_loading.bolo_data_loading as bolo_data_loading

class Bolo:

    def __init__(self, bolo_name):
        self.bolo_params = importlib.import_module("simulation.timestream_simulation.bolo_params." + bolo_name).bolo_params
        #Getting the beam profile and the del_beta
        self.beam = self.get_real_beam()
        self.axis_spin, self.axis_prec, self.axis_rev = self.get_initial_axes()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for a given bolo with any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    #@profile
    def simulate_timestream_beamed(self, segment, sky_map, out_dir):

        t_start = scan_params.t_segment*segment

        nsamples = int(1000.0*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate

        npix = hp.nside2npix(scan_params.nside)

        pad = self.del_beta.size/2
        nsamples = int(1000.0*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate + 2*pad

        npix = hp.nside2npix(scan_params.nside)

        rot_qt = self.generate_quaternion(t_start)

        pol_ang = self.get_pol_ang(rot_qt, None) 

        matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
        matrix.data.value[:, 0, 0, 0] = 0.5
        matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol_ang)
        matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol_ang)

        P = ProjectionOperator(matrix)

        signal = np.zeros(nsamples - 2*pad)

        time_start = time.time()

        for i in range(self.del_beta.size):
            print self.del_beta[i]
            sys.stdout.flush()
            v_init = self.get_initial_vec(self.del_beta[i])
            v = quaternion.transform(rot_qt, v_init)
            if scan_params.gal_coords:
                rot = hp.Rotator(coord=['E', 'G'])
                theta, phi = hp.vec2ang(v)
                theta_gal, phi_gal = rot(theta, phi)
                v = hp.ang2vec(theta_gal, phi_gal)
                del theta_gal
                del phi_gal
            hit_pix = hp.vec2pix(scan_params.nside, v[...,0], v[...,1], v[...,2])
            if self.del_beta[i] == 0.0:
                np.save(os.path.join(out_dir, "vector"), v[pad:-pad][::scan_params.oversampling_rate])
                hitpix_central = hit_pix
            del v
            P.matrix.data.index[:, 0] = hit_pix
            if scan_params.do_filtering:
                signal_int = np.convolve(P(sky_map.T), self.beam_kernel[i], mode = 'valid')
                signal += filters.filter_butter(signal_int, 1000.0/scan_params.t_sampling, scan_params.filter_cutoff)
            else:
                signal += np.convolve(P(sky_map.T), self.beam_kernel[i], mode = 'valid')

        beam_sum = np.sum(self.beam_kernel)
        signal /= beam_sum

        time_end = time.time()

        del hit_pix

        print "Average time taken per beam slice :", (time_end - time_start)/self.del_beta.size, "seconds"
        sys.stdout.flush()

        if scan_params.add_noise:
            signal += np.random.normal(scale=scan_params.noise_level, size=signal.size)

        np.save(os.path.join(out_dir, "pol_ang"), pol_ang[pad:-pad][::scan_params.oversampling_rate])
        np.save(os.path.join(out_dir, "ts_signal"), signal[::scan_params.oversampling_rate])
        del pol_ang
        del signal

        hitmap = self.get_hitmap(hitpix_central)

        return hitmap 

    def get_real_beam(self):
        beam = hp.mrdfits(self.bolo_params.beam_file)
        return beam

    #@profile
    def get_hitmap(self, hit_pix):
        nsamples = hit_pix.size
        npix = hp.nside2npix(scan_params.nside)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
        matrix.data.value = 1
        P = ProjectionOperator(matrix, shapein=npix, shapeout=nsamples)
        matrix.data.index = hit_pix[..., None]
        P.matrix = matrix
        hitmap = P.T(np.ones(nsamples, dtype=np.float32))

        return hitmap
    
    #@profile
    def generate_quaternion(self, t_start):

        if beam_params.do_pencil_beam:
            nsamples = int(1000.0*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate
            pad = 0
        else:
            pad = self.del_beta.size/2
            nsamples = int(1000.0*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate + 2*pad

        t_steps = t_start + 0.001*(scan_params.t_sampling/scan_params.oversampling_rate)*np.arange(-pad, nsamples - pad)

        w_spin = -2*np.pi/scan_params.t_spin
        w_prec = -2*np.pi/scan_params.t_prec
        w_rev = -2*np.pi/scan_params.t_year

        r_total = quaternion.multiply(quaternion.make_quaternion(w_rev*t_steps, self.axis_rev), quaternion.multiply(quaternion.make_quaternion(w_prec*t_steps, self.axis_prec), quaternion.make_quaternion(w_spin*t_steps, self.axis_spin)))

        return r_total

    #@profile
    def get_pol_ang(self, rot_qt, v_dir=None):

        pad = self.del_beta.size/2
        n_steps = int(1000.0*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate + 2*pad

        pol_init = np.deg2rad(self.bolo_params.pol_ang)
        x_axis = np.array([0.0, 1.0, 0.0])

        pol_vec = quaternion.transform(rot_qt, np.tile(x_axis, n_steps).reshape(-1,3))

        if v_dir is None:
            v_init = self.get_initial_vec(0.0)
            v_dir = quaternion.transform(rot_qt, v_init)

        theta, phi = hp.vec2ang(v_dir)

        x_local = np.array(zip(np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)))
        y_local = np.array(zip(-np.sin(phi), np.cos(phi), np.zeros(phi.size)))

        proj_x = np.sum(pol_vec*x_local, axis=-1)
        proj_y = np.sum(pol_vec*y_local, axis=-1)

        pol_ang = (np.arctan2(proj_y, proj_x) + pol_init) % np.pi 
        #pol_ang *= 2.0

        return pol_ang 

    def get_initial_vec(self, del_beta):
        alpha = np.deg2rad(scan_params.alpha)
        beta = np.deg2rad(scan_params.beta)
        del_x = np.deg2rad(self.bolo_params.del_x/60.0)
        del_y = np.deg2rad(self.bolo_params.del_y/60.0)
        del_beta_rad = np.deg2rad(del_beta/60.0)

        u_view = np.array([np.cos(alpha + beta + del_beta_rad), 0.0, np.sin(alpha + beta + del_beta_rad)])

        return u_view

    def get_initial_axes(self):
        alpha = np.deg2rad(scan_params.alpha)
        beta = np.deg2rad(scan_params.beta)

        axis_spin = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        axis_prec = np.array([1.0, 0.0, 0.0])
        axis_rev = np.array([0.0, 0.0, 1.0])

        return axis_spin, axis_prec, axis_rev


def calculate_params():

    scan_params.t_sampling = 1000.0/scan_params.sampling_rate

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

def display_params():
    print "alpha :", scan_params.alpha, " degrees"
    print "beta :", scan_params.beta, " degrees"
    print "Mode :", scan_params.mode
    t_flight = scan_params.t_segment*len(scan_params.segment_list)
    print "T flight :", t_flight/60.0/60.0, "hours /", t_flight/60.0/60.0/24, "days"
    print "T segment :", scan_params.t_segment/60.0/60.0, "hours /", scan_params.t_segment/60.0/60.0/24, "days"
    print "T precession :", scan_params.t_prec/60.0/60.0, "hours"
    print "T spin :", scan_params.t_spin, " seconds"
    print "T sampling :", scan_params.t_sampling, " milli-seconds"
    print "Scan frequency :", 1000.0/scan_params.t_sampling, "Hz"
    print "Theta co :", scan_params.theta_co, " arcmin"
    print "Theta cross :", scan_params.theta_cross, " arcmin"
    print "Scan resolution for beam integration :", scan_params.scan_resolution, "arcmin"
    print "Beam resolution :", beam_params.beam_resolution, "arcmin"
    print "Pixel size for NSIDE =", scan_params.nside, ":", hp.nside2resol(scan_params.nside, arcmin=True), "arcmin"
    print "Filter cutoff :", scan_params.filter_cutoff
    n_steps = int(1000*scan_params.t_segment/scan_params.t_sampling)*scan_params.oversampling_rate
    print "#Samples per segment :", n_steps
    print "Size of signal array :", n_steps*8.0/1024/1024, "MB"
    print "Estimated use of memory :", 15*n_steps*8.0/1024/1024, "MB"
    print "# of processes required :", len(scan_params.bolo_names)*len(scan_params.segment_list)

#@profile
def get_scanned_map(sky_map, hitmap):
    valid = hitmap>0
    sky_map[...,~valid] = np.nan
    return sky_map

def get_sky_map(bolo_name):
    sky_map = np.array(hp.read_map(scan_params.input_maps[bolo_name], field=(0,1,2)))
    nside = hp.get_nside(sky_map)
    if scan_params.nside != nside:
        print "NSIDE of map does not match with given NSIDE. Running with map NSIDE"
        scan_params.nside = nside
    if scan_params.do_only_T:
        sky_map[1:,...] = 0.0
    return sky_map

def make_output_dirs(out_dir, scan_params):
    os.makedirs(out_dir)
    param_dir = os.path.join(out_dir, "params")
    bolo_param_dir = os.path.join(param_dir, "bolo_params")
    os.makedirs(param_dir)
    os.makedirs(bolo_param_dir)
    shutil.copy("default_params.py", param_dir)
    shutil.copy("custom_params.py", param_dir)
    shutil.copy("sim_timestream_pol.py", out_dir)
    for bolo in scan_params.bolo_names:
        os.makedirs(os.path.join(out_dir, bolo))
        shutil.copy("bolo_params/"+bolo+".py", bolo_param_dir)
        for segment in scan_params.segment_list:
            segment_name = str(segment+1).zfill(4)
            os.makedirs(os.path.join(out_dir, bolo, segment_name))
            

def run_serial():

    hitmap = np.zeros(hp.nside2npix(scan_params.nside))

    time_stamp = get_time_stamp()
    out_dir = os.path.join(scan_params.global_output_dir, "scanning", time_stamp)
    make_output_dirs(out_dir, scan_params)
    
    calculate_params()

    display_params()

    for bolo_name in scan_params.bolo_names:
        sky_map = get_sky_map(bolo_name) 
        bolo = Bolo(bolo_name)
        print "Doing Bolo : ", bolo_name
        for segment in scan_params.segment_list:
            out_dir_local = os.path.join(out_dir, bolo_name, str(segment+1).zfill(4))
            segment_name = str(segment+1).zfill(4)
            print "Segment : ", segment_name
            hitmap += bolo.simulate_timestream(segment, sky_map, out_dir_local)



    scanned_map = get_scanned_map(sky_map, hitmap)
    hp.write_map(os.path.join(out_dir, "scanned_map.fits"), scanned_map)
    hp.write_map(os.path.join(out_dir, "hitmap_in.fits"), hitmap)

    return hitmap 


#@profile
def run_mpi():
    time_start = time.time()   

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    time_stamp = [None]
    if rank is 0:
        time_stamp[0] = get_time_stamp()
    time_stamp = comm.bcast(time_stamp, root=0)

    out_dir = os.path.join(scan_params.global_output_dir, "scanning", time_stamp[0])

    calculate_params()

    if rank is 0:
        display_params()
        make_output_dirs(out_dir, scan_params)

    comm.Barrier()

    hitmap_local = np.zeros(hp.nside2npix(scan_params.nside), dtype=np.float32)

    local_bolo_list, local_segment_list = bolo_data_loading.get_local_bolo_segment_list(rank, size, scan_params)

    prev_bolo_name = None 

    time_segment_start = time.time()

    for bolo_name, segment in zip(local_bolo_list, local_segment_list):
        if bolo_name != prev_bolo_name:
            sky_map = get_sky_map(bolo_name)
            bolo = Bolo(bolo_name)
        print "Doing Bolo :", bolo_name, " Segment :", segment+1, " Rank :", rank
        sys.stdout.flush()
        out_dir_local = os.path.join(out_dir, bolo_name, str(segment+1).zfill(4))
        if beam_params.do_pencil_beam:
            hitmap_local += bolo.simulate_timestream(segment, sky_map, out_dir_local)
        else:
            hitmap_local += bolo.simulate_timestream_beamed(segment, sky_map, out_dir_local)

    time_segment_stop = time.time()
    print "Average time taken per segment :", (time_segment_stop - time_segment_start)/len(local_bolo_list), "s"
    sys.stdout.flush()

    if rank is 0:
        hitmap = np.zeros(hp.nside2npix(scan_params.nside), dtype=np.float32)
    else:
        hitmap = None

    comm.Reduce(hitmap_local, hitmap, MPI.SUM, 0)

    if rank is 0:
        valid = hitmap>0
        sky_map[...,~valid] = np.nan
        hp.write_map(os.path.join(out_dir, "scanned_map.fits"), sky_map)
        hp.write_map(os.path.join(out_dir, "hitmap_in.fits"), hitmap)

    time_end = time.time() 
    if rank is 0:
        print "Total time taken :", (time_end - time_start), "seconds"
        print "Scan time stamp :", time_stamp[0]


if __name__=="__main__":
    from custom_params import scan_params, beam_params
    action = sys.argv[1]

    if action=='display_params':
        calculate_params()
        display_params()

    if action=='run_mpi':
        from mpi4py import MPI
        run_mpi()

    if action=='run_serial':
        run_serial()
