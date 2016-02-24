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

def estimate_power_spectrum(sky_map_masked, hitmap, lmax=2000):
    spectrum = hp.anafast((sky_map_masked[0].filled(), sky_map_masked[1].filled(), sky_map_masked[2].filled()), lmax=lmax)
    f_sky = float(np.sum(hitmap>0))/hitmap.size
    return np.array(spectrum)[...,2:]/f_sky

map_file = "/global/homes/b/banerji/simulation/output/reconstructing/2016_02_18__13_38_03/reconstructed_map.fits"
hitmap_file = "/global/homes/b/banerji/simulation/output/reconstructing/2016_02_18__13_38_03/hitmap_out.fits"
spectra_file = "/global/homes/b/banerji/simulation/output/reconstructing/2016_02_18__13_38_03/power_spectra"

sky = hp.read_map(map_file, field=(0,1,2))
hitmap = hp.read_map(hitmap_file)

sky_masked = mask_map(sky, hitmap)

spectra = estimate_power_spectrum(sky_masked, hitmap)

np.save(spectra_file, spectra)
