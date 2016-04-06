#!/usr/bin/env python

import numpy as np
import healpy as hp
import simulation.lib.numericals.filters as fl
from simulation.params.custom_params import global_paths
import os

def mask_map(sky_map, binary_mask=None, pol=False):
    if binary_mask is None:
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

def estimate_cl(sky_map, lmax, binary_mask=None, fwhm=0.0, pol=False):

    sky_map_masked = mask_map(sky_map, binary_mask=binary_mask, pol=pol)

    if pol:
        binary_mask = np.logical_not(sky_map_masked[0].mask)
    else:
        binary_mask = np.logical_not(sky_map_masked.mask)
    f_sky = float(np.sum(binary_mask))/binary_mask.size

    Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=pol)

    if pol:
        spectra = hp.anafast((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
        spectra /= f_sky*Bl.T**2
    else:
        spectra = hp.anafast(sky_map_masked.filled(), lmax=lmax)
        spectra /= f_sky*Bl**2

    spectra = np.array(spectra)
    return spectra

def wiener_filter_for_alm(alm, lmax, fwhm=0.0, f_sky=1.0, sky_prior=None):
    if sky_prior is None:
        spectra_th = np.load("/global/homes/b/banerji/simulation/spectra/r_001/unlensed_cls.npy")[0,:lmax+1]
    else:
        spectra_th = estimate_cl(sky_prior, lmax, fwhm=fwhm)

    Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=False)
    spectra_ob = hp.alm2cl(alm, lmax_out=lmax)
    spectra_ob /= f_sky*Bl**2

    filter_response = spectra_ob/spectra_th
    filter_response[:2] = 1.0
    filter_response = np.sqrt(filter_response)

    return filter_response

"""
def get_wiener_filter(sky_map, lmax, alm=None, binary_mask=None, fwhm=0.0):
    spectra_th = np.load("/global/homes/b/banerji/simulation/spectra/r_001/unlensed_cls.npy")[0,:lmax+1]

    if sky_map is None:
        Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=False)
        f_sky = float(np.sum(binary_mask))/binary_mask.size
        spectra_ob = hp.alm2cl(alm, lmax_out=lmax)
        spectra_ob /= f_sky*Bl**2
    else:
        spectra_ob = estimate_cl(sky_map, lmax, fwhm=fwhm, binary_mask=binary_mask, pol=False) 

    filter_response = spectra_ob/spectra_th
    filter_response[:2] = 1.0

    return filter_response
"""

def estimate_alm(sky_map, lmax, binary_mask=None, pol=False):
    sky_map_masked = mask_map(sky_map, binary_mask, pol=pol)
    if pol:
        alm = hp.map2alm((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
    else:
        alm = hp.map2alm(sky_map_masked.filled(), lmax=lmax)

    return alm

def deconvolve_alm(alms, lmax=None, fwhm=0.0, f_sky=1.0, pol=False, wiener=True, sky_prior=None):
    if fwhm == 0.0:
        return alms

    retalm = []

    factor = (2.0*np.sqrt(2.0*np.log(2.0)))
    sigma = fwhm/factor
    
    if lmax is None:
        lmax = hp.Alm.getlmax(len(alms[0] if pol else alms), None)
    
    if wiener:
        wiener_filter = wiener_filter_for_alm(alms[0] if pol else alms, lmax, fwhm=fwhm, f_sky=f_sky, sky_prior=sky_prior)
        wiener_smooth=wiener_filter
        wiener_smooth = fl.filter_butter(wiener_filter, lmax, 100)
        wiener_smooth[:1500] = 1.0

    if pol:
        for ialm, alm in enumerate(alms):
            ell = np.arange(lmax + 1)
            if ialm >= 1:
                s = 2
            else:
                s = 0
            fact = np.exp(0.5*(ell*(ell + 1) - s**2)*sigma**2)
            if wiener:
                fact /= wiener_smooth
            res = hp.almxfl(alm, fact, inplace=False)
            retalm.append(res)
        retalm = np.array(retalm)
    else:
        lmax = hp.Alm.getlmax(len(alms), None)
        ell = np.arange(lmax + 1)
        fact = np.exp(0.5*(ell*(ell + 1))*sigma**2)
        if wiener:
            fact /= wiener_smooth
        retalm = hp.almxfl(alms, fact, inplace=False)

    return retalm 


def deconvolve_map(map_in, fwhm, lmax=None, binary_mask=None, pol=False, wiener=True, sky_prior=None):
    #if fwhm == 0:
    #    return map_in

    if binary_mask is None:
        binary_mask = np.logical_not(np.isnan(map_in))
    f_sky = float(np.sum(binary_mask))/binary_mask.size
    alm_in = estimate_alm(map_in, lmax, binary_mask, pol)
    alm_dec = deconvolve_alm(alm_in, fwhm=fwhm, f_sky=f_sky, pol=pol, wiener=True, sky_prior=sky_prior)
    map_dec = hp.alm2map(alm_dec, nside=hp.get_nside(map_in))
    return map_dec
