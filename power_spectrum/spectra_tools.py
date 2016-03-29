#!/usr/bin/env python

import numpy as np
import healpy as hp
from simulation.params.custom_params import global_paths
import os

def mask_map(sky_map, binary_mask=None):
    if binary_mask==None:
        binary_mask = np.logical_not(np.isnan(sky_map[0]))

    sky_map_masked = hp.ma(sky_map)

    sky_map_masked[0].mask = np.logical_not(binary_mask)
    sky_map_masked[1].mask = np.logical_not(binary_mask)
    sky_map_masked[2].mask = np.logical_not(binary_mask)

    return sky_map_masked

def estimate_cl(sky_map, lmax=2000, fwhm=8.0):
    binary_mask = np.logical_not(np.isnan(sky_map[0]))
    sky_map_masked = mask_map(sky_map, binary_mask)
    spectra = hp.anafast((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
    f_sky = float(np.sum(binary_mask))/binary_mask.size
    ell = np.arange(spectra[0].size)
    factor = np.sqrt(8*np.log(2))
    sigma = np.radians(fwhm/60.0/factor)
    Bl = np.exp(-0.5*ell*(ell+1)*sigma**2)
    spectra = np.array(spectra)
    spectra /= f_sky*Bl**2
    return spectra[...,2:]

def estimate_alm(sky_map, lmax=2000):
    binary_mask = np.logical_not(np.isnan(sky_map[0]))
    sky_map_masked = mask_map(sky_map, binary_mask)
    alm = hp.map2alm((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)

def plot_theoretical_bb(r=['01', '001', '0001'], lmax=2000):
    ell = np.arange(2, lmax+1)
    for r_value in r:
        spectra = np.load("../spectra/" + r_value + "/unlensed_cls.py")[2, 2:lmax+1]
        loglog(ell, ell*(ell+1)*spectra/2/np.pi)

def deconvolve_alm(alms, fwhm, pol=True):
    if fwhm == 0:
        return alm

    retalm = []
    sigma = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
    
    if pol:
        for ialm, alm in enumerate(alms):
            lmax = hp.Alm.getlmax(len(alm), None)
            ell = np.arange(lmax + 1)
            if ialm >= 1:
                s = 2
            else:
                s = 0
            fact = np.exp(0.5*(ell*(ell + 1) - s**2)*sigma**2)
            res = hp.almxfl(alm, fact)
            retalm.append(res)
    else:
        lmax = hp.Alm.getlmax(len(alms), None)
        ell = np.arange(lmax + 1)
        fact = np.exp(0.5*(ell*(ell + 1))*sigma**2)
        return hp.almxfl(alms, fact)

    return retalm
