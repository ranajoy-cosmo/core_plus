#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys, importlib
import pyoperators as po
from pysimulators import BeamGaussian
from simulation.lib.plotting.my_imshow import new_imshow
from custom_settings import settings
from simulation.bolo.custom_settings import pointing_params
bolo_params = importlib.import_module("simulation.bolo.bolo_params.0001").bolo_params


def generate_pointing(del_beta):
    beta = np.radians(bolo_params.beta + del_beta/60.0)
    u_init = np.array([np.cos(beta), 0.0, np.sin(beta)])
    w_prec = 2*np.pi/pointing_params.t_prec
    w_spin = 2*np.pi/pointing_params.t_spin
    t_steps = 0.001*(pointing_params.t_sampling/pointing_params.oversampling_rate)*np.arange(-n, n+1)
    if del_beta==0.0:
        print "Initial vector :", u_init
        print "t steps :", t_steps
    R = po.Rotation3dOperator("XY'X''", -1.0*w_prec*t_steps, -1.0*np.full(2*n+1, np.radians(bolo_params.alpha)), -w_spin*t_steps)
    v = R*u_init
    return v

settings.nside = 4096*4

beam = BeamGaussian(np.radians(settings.fwhm_major/60.0))
beam_healpix = beam.healpix(settings.nside)

print "Beta :", bolo_params.beta
print "Alpha :", bolo_params.alpha
print "FWHM :", settings.fwhm_major
print "NSIDE :", settings.nside

pointing_healpix = np.zeros(beam_healpix.size)

if pointing_params.mode==1:
    settings.t_spin = 360.0*60.0*np.sin(bolo_params.beta)*settings.t_sampling/1000.0/settings.theta_co
    settings.t_prec = 360.0*60.0*np.sin(bolo_params.alpha)*settings.t_spin/settings.theta_cross

if pointing_params.mode==2:
    pointing_params.theta_cross = 360.0*60.0*np.sin(bolo_params.alpha)*pointing_params.t_spin/pointing_params.t_prec
    pointing_params.theta_co = 360*60*np.sin(bolo_params.beta)*pointing_params.t_sampling/1000.0/pointing_params.t_spin

del_beta = pointing_params.theta_co/pointing_params.oversampling_rate
print "T Precession :", pointing_params.t_prec/60.0/60.0
print "T Spin :", pointing_params.t_spin
print "Theta-co :", pointing_params.theta_co
print "Theta-cross :", pointing_params.theta_cross
print "Oversampling rate :", pointing_params.oversampling_rate
print "Del beta :", del_beta

n = int(2*pointing_params.oversampling_rate*settings.fwhm_major/pointing_params.theta_co)
#n = int(1000*pointing_params.t_spin/pointing_params.t_sampling)
del_betas = np.arange(-n, n+1)*del_beta
print del_betas

beam_kernel = np.zeros(2*n+1)

for del_beta in del_betas:
    v = generate_pointing(del_beta)
    lat, lon = hp.vec2ang(v)
    pix = hp.vec2pix(settings.nside, v[...,0], v[...,1], v[...,2])
    
    beam_row = hp.get_interp_val(beam_healpix, lat, lon)
    beam_kernel = np.vstack((beam_kernel, beam_row))
    
    """
    print pix
    if del_beta==0.0:
        #v = generate_pointing(del_beta)
        #pix = hp.vec2pix(settings.nside, v[...,0], v[...,1], v[...,2])
        lat, lon = hp.pix2ang(settings.nside, pix)
        print "Lat :", lat
        print "Lon :", lon
    pointing_healpix[pix] = 1
    """

beam_kernel = beam_kernel[1:]
print beam_kernel.shape


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

plot_beam()
