#!/usr/bin/env python

import numpy as np
import healpy as hp
from simulation.params.custom_params import global_paths
import os

def mask_map(sky_map, binary_mask=None, pol=True):
    if binary_mask==None:
        if pol:
            binary_mask = np.logical_not(np.isnan(sky_map[0]))
        else:
            binary_mask = np.logical_not(np.isnan(sky_map))

    sky_map_masked = hp.ma(sky_map)

    if pol:
        sky_map_masked[0].mask = np.logical_not(binary_mask)
        sky_map_masked[1].mask = np.logical_not(binary_mask)
        sky_map_masked[2].mask = np.logical_not(binary_mask)
    else:
        sky_map_masked.mask = np.logical_not(binary_mask)

    return sky_map_masked

def estimate_cl(sy_map, fwhm=0.0, binary_mask=None, lmax=5000, pol=True):

    sky_map_masked = mask_map(sky_map, binary_mask, pol=pol)

    if pol:
        binary_mask = np.logical_not(sky_map_masked[0].mask)
    else:
        binary_mask = np.logical_not(sky_map_masked.mask)
    f_sky = float(np.sum(binary_mask))/binary_mask.size
    #ell = np.arange(lmax+1)
    #factor = np.sqrt(8*np.log(2))
    #sigma = fwhm/factor
    Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=pol)

    if pol:
        spectra = hp.anafast((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
    else:
        spectra = hp.anafast(sky_map_masked.filled(), lmax=lmax)
        #Bl = np.exp(-0.5*ell*(ell+1)*sigma**2)

    spectra = np.array(spectra)
    spectra /= f_sky*Bl**2
    return spectra

def estimate_alm(sky_map, binary_mask=None, lmax=5000, pol=True):
    sky_map_masked = mask_map(sky_map, binary_mask, pol=pol)
    if pol:
        alm = hp.map2alm((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
    else:
        alm = hp.map2alm(sky_map_masked.filled(), lmax=lmax)

    return alm

def plot_theoretical_bb(r=['01', '001', '0001'], lmax=2000):
    ell = np.arange(2, lmax+1)
    for r_value in r:
        spectra = np.load("../spectra/" + r_value + "/unlensed_cls.py")[2, 2:lmax+1]
        loglog(ell, ell*(ell+1)*spectra/2/np.pi)

def deconvolve_alm(alms, fwhm=0.0, pol=True):
    if fwhm == 0.0:
        return alms

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
            res = hp.almxfl(alm, fact, inplace=False)
            retalm.append(res)
    else:
        lmax = hp.Alm.getlmax(len(alms), None)
        ell = np.arange(lmax + 1)
        fact = np.exp(0.5*(ell*(ell + 1))*sigma**2)
        return hp.almxfl(alms, fact, inplace=False)

def get_weiner_filter(sky_map, fwhm=0.0, binary_map=None, lmax=5000):
    spectra_th = np.load("/global/homes/b/banerji/simulation/spectra/r_001/unlensed_cls.npy")[0,:lmax+1]

    spectra_ob = estimate_cl(sky_map, fwhm, binary_map, lmax, False) 

    print spectra_th.shape
    print spectra_ob.shape

    filter_response = spectra_ob/spectra_th

    return filter_response


def deconvolve_map(map_in, fwhm, binary_mask=None, lmax=5000, pol=True):
    #if fwhm == 0:
    #    return map_in

    alm_in = estimate_alm(map_in, binary_mask, lmax, pol)
    alm_dec = deconvolve_alm(alm_in, fwhm=fwhm, pol=pol)
    map_dec = hp.alm2map(alm_dec, nside=hp.get_nside(map_in))
    return map_dec
