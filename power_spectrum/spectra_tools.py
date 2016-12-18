#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import simulation.lib.numericals.filters as fl
from simulation.global_config import global_paths
import os

def mask_map(sky_map, binary_mask=None, pol=True, ret_mask=False, fill_zeros=False):

    if binary_mask is None:
        binary_mask = np.logical_not(np.isnan(sky_map[0] if pol else sky_map))

    if fill_zeros:
        if pol:
            sky_map[0][np.isnan(sky_map[0])] = 0.0
            sky_map[1][np.isnan(sky_map[1])] = 0.0
            sky_map[2][np.isnan(sky_map[2])] = 0.0
        else:
            sky_map[np.isnan(sky_map)] = 0.0

    sky_map_masked = hp.ma(sky_map)

    if pol:
        sky_map_masked[0].mask = np.logical_not(binary_mask)
        sky_map_masked[1].mask = np.logical_not(binary_mask)
        sky_map_masked[2].mask = np.logical_not(binary_mask)
    else:
        sky_map_masked.mask = np.logical_not(binary_mask)

    if ret_mask:
        return sky_map_masked, binary_mask
    else:
        return sky_map_masked


def mask_map_apodised(sky_map, apodised_mask, pol=True):
    nside = hp.get_nside(sky_map)
    if pol:
        masked_map = np.empty((3, 12*nside**2))
        for i in range(3):
            masked_map[i] = sky_map[i]*apodised_mask
    else:
        masked_map = sky_map*apodised_mask

    f_sky = np.mean(apodised_mask**2)

    return masked_map, f_sky


def estimate_cl(sky_map, lmax, binary_mask=None, fwhm=0.0, pol=True):

    sky_map_masked, binary_mask = mask_map(sky_map, binary_mask=binary_mask, pol=pol, ret_mask=True, fill_zeros=True)

    f_sky = float(np.sum(binary_mask))/binary_mask.size

    Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=pol)

    if pol:
        spectra = hp.anafast((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)[:4]
        spectra = np.array(spectra)
        spectra /= f_sky*Bl.T**2
    else:
        spectra = hp.anafast(sky_map_masked.filled(), lmax=lmax)
        spectra /= f_sky*Bl**2

    spectra = np.array(spectra)
    return spectra


def estimate_alm(sky_map, lmax, binary_mask=None, pol=False):
    sky_map_masked = mask_map(sky_map, binary_mask, pol=pol)
    if pol:
        alm = hp.map2alm((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
    else:
        alm = hp.map2alm(sky_map_masked.filled(), lmax=lmax)

    return alm


def wiener_filter_for_alm(alm, lmax=None, fwhm=0.0, f_sky=1.0, sky_prior=None):

    if lmax is None:
        lmax = hp.Alm.getlmax(len(alm), None)

    if sky_prior is None:
        spectra_th = np.load("/global/homes/b/banerji/simulation/spectra/r_001/lensedtot_cls.npy")[0,:lmax+1]
    else:
        spectra_th = estimate_cl(sky_prior, lmax, fwhm=fwhm, pol=False)

    Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=False)
    spectra_ob = hp.alm2cl(alm, lmax_out=lmax)
    spectra_ob /= f_sky*Bl**2

    filter_response = spectra_ob/spectra_th
    filter_response[:2] = 1.0
    filter_response = np.sqrt(filter_response)

    return filter_response


def deconvolve_alm(alms, lmax=None, fwhm_in=0.0, fwhm_out=0.0, f_sky=1.0, pol=False, wiener=True, sky_prior=None):

    fwhm = np.sqrt(fwhm_in**2 - fwhm_out**2)
    if fwhm == 0.0:
        return alms
    factor = (2.0*np.sqrt(2.0*np.log(2.0)))
    sigma = fwhm/factor
    
    retalm = []

    if lmax is None:
        lmax = hp.Alm.getlmax(len(alms[0] if pol else alms), None)
    
    if wiener:
        wiener_filter = wiener_filter_for_alm(alms[0] if pol else alms, lmax, f_sky=f_sky, sky_prior=sky_prior, fwhm=fwhm_in)
        wiener_smooth = fl.filter_butter(wiener_filter, lmax, 100)
        #wiener_smooth[:1500] = 1.0

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
    else:
        lmax = hp.Alm.getlmax(len(alms), None)
        ell = np.arange(lmax + 1)
        fact = np.exp(0.5*(ell*(ell + 1))*sigma**2)
        if wiener:
            fact /= wiener_smooth
        retalm = hp.almxfl(alms, fact, inplace=False)

    return retalm 


def deconvolve_map(map_in, fwhm_in=0.0, fwhm_out=0.0, lmax=None, binary_mask=None, pol=False, wiener=True, sky_prior=None):
    
    if fwhm_in == fwhm_out:
        return map_in

    if lmax is None:
        lmax = 3*hp.get_nside(map_in) - 1

    if binary_mask is None:
        binary_mask = np.logical_not(np.isnan(map_in[0] if pol else map_in))
    f_sky = float(np.sum(binary_mask))/binary_mask.size

    alm_in = estimate_alm(map_in, lmax, binary_mask, pol)
    alm_dec = deconvolve_alm(alm_in, fwhm_in=fwhm_in, fwhm_out=fwhm_out, f_sky=f_sky, pol=pol, wiener=True, sky_prior=sky_prior)
    map_dec = hp.alm2map(alm_dec, nside=hp.get_nside(map_in), pol=pol)
    
    return map_dec

def get_central_mask(nside, radius, deg=True):

    pix = np.arange(12*nside**2)
    theta_c, phi_c = np.pi/2, 0.0

    pix_dist = hp.rotator.angdist([theta_c, phi_c], hp.pix2ang(nside, pix))

    if deg:
        radius = np.radians(radius)

    mask = pix_dist<radius

    return mask


def get_temp_gradient_map(T_sky, fwhm_in=0.0, fwhm_out=0.0, nside=None, lmax=None):
    if nside==None:
        nside = hp.get_nside(T_sky)
    if lmax==None:
        lmax = 3*nside - 1

    alm = estimate_alm(T_sky, lmax)
    alm_dec = deconvolve_alm(alm, lmax, fwhm_in, fwhm_out)
    grad_sky_T = hp.alm2map_der1(alm_dec, nside, lmax)

    return grad_sky_T[1:]


def bin_spectra_old(spectra, bins, pol=False):
    n_bins = len(bins) - 1
    ell = np.arange(bins[-1] + 1)
    ell_bins = np.empty(n_bins)
    if pol:
        binned_spectra = np.empty((3, n_bins))
    else:
        binned_spectra = np.empty(n_bins)

    for i in range(n_bins):
        bin_low = bins[i]
        bin_high = bins[i+1]
        if pol:
            for j in range(3):
                binned_spectra[j, i] = np.mean(spectra[j, bin_low:bin_high])
        else:
            binned_spectra[i] = np.mean(spectra[bin_low:bin_high])
        ell_bins[i] = np.mean(ell[bin_low:bin_high])

    return ell_bins, binned_spectra

def bin_spectra(spectra, lmax, pol=False):
    if pol:
        binned_spectra = np.empty((3, lmax+1))
    else:
        binned_spectra = np.empty(lmax+1)

    for ell in range(lmax+1):
        ell_low = ell - 0.1*ell
        ell_high = ell + 0.1*ell
        if pol:
            for j in range(3):
                binned_spectra[j, ell] = np.mean(spectra[j, ell_low:ell_high])
        else:
            binned_spectra[ell] = np.mean(spectra[ell_low:ell_high])

    return binned_spectra

def bin_error_bar(error, bins, pol=False):
    n_bins = len(bins) - 1
    ell = np.arange(bins[-1] + 1)
    bin_width = np.empty(n_bins)
    if pol:
        binned_error = np.empty((3, n_bins))
    else:
        binned_error = np.empty(n_bins)

    for i in range(n_bins - 1):
        bin_low = bins[i]
        bin_high = bin[i+1]
        bin_width[i] = bin_high - bin_low
        if pol:
            for j in range(3):
                binned_error[j, i] = np.mean(error[j, bin_low:bin_high])/np.sqrt(bin_width[i])
        else:
            binned_error[i] = np.mean(error[bin_low:bin_high])/np.sqrt(bin_width[i])

    return bin_width, binned_error

def TQU_TO_TEB_maps(sky_TQU, lmax):
    nside = hp.get_nside(sky_TQU)

    alm_in = hp.map2alm(sky_TQU, lmax, pol=True)

    sky_TEB = np.empty((3, 12*nside**2))
    for i in range(3):
        sky_TEB[i] = hp.alm2map(alm_in[i], nside, pol=False)

    return sky_TEB

def TEB_to_TQU_maps(sky_TEB, lmax):
    nside = hp.get_nside(sky_TEB)

    alm = np.empty((3, (lmax + 1)**2))
    
    for i in range(3):
        alm[i] = hp.map2alm(sky_TEB[i], lmax, pol=False)

    sky_TQU = hp.alm2map(alm, nside, lmax, pol=True)

    return sky_TQU
