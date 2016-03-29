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

def deconvolve_alm(alm, fwhm):
    if fwhm == 0:
        return alm

    retalm = []
    sigma = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
    
    for ialm, alm in enumerate(alm):
        lmax = hp.Alm.getlmax(len(alm), None)
        ell = np.arange(lmax + 1)
        if ialm >= 1:
            s = 2
        else:
            s = 0
        fact = np.exp(0.5*(ell*(ell + 1) - s**2)*sigma**2)
        res = hp.almxfl(alm, fact)
        retalm.append(res)

    return retalm

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()

#dir = global_paths.maps_dir
#map_dirs = ["2016_03_16__00_32_24", "2016_03_16__10_30_53", "2016_03_16__11_41_40", "2016_03_16__12_28_20"]#, "2016_03_16__12_40_51"]
map_dirs = "fourth_set/2016_02_24__21_24_47"

#dir = os.path.join(global_paths.output_dir, "reconstructing", map_dirs[rank])
dir = os.path.join("/scratch1/scratchdirs/banerji/core_long_term_output/reconstruction", map_dirs)#[rank])

map_file = os.path.join(dir, "leakage_subtracted.fits")
spectra_file = os.path.join(dir, "power_spectra_sub")

sky = hp.read_map(map_file, field=(0,1,2))
binary_mask = np.logical_not(np.isnan(sky[0]))

sky_masked = mask_map(sky, binary_mask)

spectra = estimate_power_spectrum(sky_masked, binary_mask, lmax=2000, fwhm=8.0)

np.save(spectra_file, spectra)
