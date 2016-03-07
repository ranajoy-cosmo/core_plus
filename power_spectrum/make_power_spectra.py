#!/usr/bin/env python

import numpy as np
import healpy as hp
from simulation.params.custom_params import global_paths
import os

def mask_map(sky_map, binary_mask):
    sky_map_masked = hp.ma(sky_map)

    sky_map_masked[0].mask = np.logical_not(binary_mask)
    sky_map_masked[1].mask = np.logical_not(binary_mask)
    sky_map_masked[2].mask = np.logical_not(binary_mask)

    return sky_map_masked

def estimate_power_spectrum(sky_map_masked, binary_mask, lmax=2000, fwhm=8.0):
    spectra = hp.anafast((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
    f_sky = float(np.sum(binary_mask))/binary_mask.size
    ell = np.arange(spectra[0].size)
    factor = np.sqrt(8*np.log(2))
    sigma = np.radians(fwhm/60.0/factor)
    Bl = np.exp(-0.5*ell*(ell+1)*sigma**2)
    spectra = np.array(spectra)
    spectra /= f_sky*Bl**2
    return spectra[...,2:]

dir = global_paths.maps_dir
#dir = os.path.join(global_paths.output_dir, "scanning", "2016_03_06__19_17_46")

map_file = os.path.join(dir, "sky_map_4096_8_001.fits")
spectra_file = os.path.join(dir, "power_spectra_8_001")

sky = hp.read_map(map_file, field=(0,1,2))
binary_mask = np.logical_not(np.isnan(sky[0]))

sky_masked = mask_map(sky, binary_mask)

spectra = estimate_power_spectrum(sky_masked, binary_mask, lmax=2000, fwhm=8.0)

np.save(spectra_file, spectra)
