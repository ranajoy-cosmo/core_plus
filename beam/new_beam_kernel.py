#!/usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import pysimulators as ps
import sys
from simulation.lib.plotting.my_imshow import new_imshow

nside = 2048
def healpix_beam():
    fwhm_major = 8.0
    beam = ps.BeamGaussian(np.radians(fwhm_major/60.0))
    beam_healpix = beam.healpix(nside)
    return beam_healpix

def generate_pointing(del_beta=0.0):
    beta = np.radians(45.0)
    del_beta = np.radians(del_beta/60.0)
    alpha = np.radians(45.0)
    t_sampling = 5.0
    oversampling_rate = 2
    t_year = 365*24*60*60.0
    t_prec = 4*24*60*60.0
    t_spin = 60.0
    n_steps = 10.0

    u_init = np.array([np.cos(beta + del_beta), 0.0, np.sin(beta + del_beta)])
    w_orbit = 2*np.pi/t_year
    w_prec = 2*np.pi/t_prec
    w_spin = 2*np.pi/t_spin
    t_steps = 0.001*(t_sampling/oversampling_rate)*np.arange(-n_steps, n_steps)
    R = po.Rotation3dOperator("XY'X''", -1.0*w_prec*t_steps, -1.0*np.full(2*n_steps, alpha), -w_spin*t_steps)
    v = R*u_init
    return v

beam_healpix = healpix_beam()
del_beta_all = 2*np.arange(-5,6)
pointing_healpix = np.zeros(beam_healpix.size)
for del_beta in del_beta_all:
    v = generate_pointing(del_beta)
    pix = hp.vec2pix(nside, v[...,0], v[...,1], v[...,2])
    pointing_healpix[pix] = 1
