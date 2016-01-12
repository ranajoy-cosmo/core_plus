#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys, importlib
import pyoperators as po
from pysimulators import BeamGaussian
from simulation.lib.plotting.my_imshow import new_imshow
from custom_settings import settings
from simulations.bolo.custom_settings import pointing_params
bolo_params = importlib.import_module("simulation.bolo.bolo_params.0001").bolo_params


beam = BeamGaussian(np.radians(settings.fwhm_major/60.0))
print "FWHM :",settings.fwhm_major
beam_healpix = beam.healpix(settings.nside)
pointing_healpix = np.zeros(beam_healpix.size)
del_beta = pointing_params.theta_co/pointing_params.oversampling_rate
print "Theta-co :", pointing_params.theta_co
print "Del beta :", del_beta
n = int(3*settings.fwhm/del_beta)
del_betas = np.arange(-n, n+1)*del_beta
print del_betas
for del_beta in del_betas:
    v = generate_pointing(del_beta)
    pix = hp.vec2pix(nside, v[...,0], v[...,1], v[...,2])
    pointing_healpix[pix] = 1


def generate_pointing(del_beta):
    u_init = np.array([np.cos(bolo_params.beta + del_beta), 0.0, np.sin(bolo_params.beta + del_beta)])
    t_spin = 360.0*60.0*np.sin(bolo_params.beta)*pointing_params.t_sampling/1000.0/pointing_params.theta_co
    t_prec = 360.0*60.0*np.sin(bolo_params.alpha)*pointing_params.t_spin/pointing_params.theta_cross
    w_prec = 2*np.pi/t_prec
    w_spin = 2*np.pi/t_spin
    t_steps = 0.001*(pointing_params.t_sampling/pointing_params.oversampling_rate)*np.arange(-n, n+1)
    R = po.Rotation3dOperator("XY'X''", -1.0*w_prec*t_steps, -1.0*np.full(2*n+1, bolo_params.alpha), -w_spin*t_steps)
    v = R*u_init
    return v

def display_beam_settings():
    if settings.do_pencil_beam:
        print "Pencil beam"
    else:
        print "Major axis(FWHM) : ", settings.fwhm_major, "arcmins"
        print "Minor axis(FWHM) : ", settings.fwhm_minor, "arcmins"
        print "Center : ", settings.center
        print "Tilt : ", settings.tilt, " degrees"
        print "Pixel size : ", settings.beam_resolution, "arcmins" 
        print "No, of pixels per fwhm (minor-axis): ", 2*settings.fwhm_minor/settings.beam_resolution

def plot_beam():
    fig, ax = plt.subplots()
    im = new_imshow(ax, beam_kernel)
    fig.colorbar(im, ax=ax)
    plt.show()
