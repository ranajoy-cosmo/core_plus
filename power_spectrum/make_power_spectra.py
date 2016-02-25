#!/usr/bin/env python

import numpy as np
import healpy as hp

def mask_map(sky_map, hitmap):
    mask = hitmap>0
    sky_map_masked = hp.ma(sky_map)

    sky_map_masked[0].mask = np.logical_not(mask)
    sky_map_masked[1].mask = np.logical_not(mask)
    sky_map_masked[2].mask = np.logical_not(mask)

    return sky_map_masked

def estimate_power_spectrum(sky_map_masked, hitmap, lmax=2000, fwhm=8.0):
    spectra = hp.anafast((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
    f_sky = float(np.sum(hitmap>0))/hitmap.size
    ell = np.arange(spectra[0].size)
    factor = np.sqrt(8*np.log(2))
    sigma = np.radians(fwhm/60.0/factor)
    Bl = np.exp(-0.5*ell*(ell+1)*sigma**2)
    spectra = np.array(spectra)
    spectra /= f_sky*Bl**2
    return spectra[...,2:]
"""
map_file = "/global/homes/b/banerji/simulation/output/reconstructing/2016_02_24__21_45_56/reconstructed_map.fits"
hitmap_file = "/global/homes/b/banerji/simulation/output/reconstructing/2016_02_24__21_45_56/hitmap_out.fits"
spectra_file = "/global/homes/b/banerji/simulation/output/reconstructing/2016_02_24__21_45_56/power_spectra"

sky = hp.read_map(map_file, field=(0,1,2))
hitmap = hp.read_map(hitmap_file)

sky_masked = mask_map(sky, hitmap)

spectra = estimate_power_spectrum(sky_masked, hitmap)

np.save(spectra_file, spectra)
"""
