#!/usr/bin/env python 

import numpy as np
import healpy as hp

LMAX = 1500
NSIDE = 1024

"""
mask = hp.read_map("hit_map.fits")
mask[mask != 0] = 1
hp.write_map("binary_mask.fits", mask)

#sky_map = hp.ud_grade(hp.read_map("/global/homes/b/banerji/leap/apps/simulated_timestreams/bolo/maps_and_spectra/sky_map_1024.fits", field = (0,1,2)), nside_out = NSIDE)
sky_map = hp.read_map("/global/homes/b/banerji/leap/apps/simulated_timestreams/bolo/maps_and_spectra/sky_map_1024.fits", field = (0,1,2))
for i in (0,1,2):
	sky_map[i][np.logical_not(mask)] = np.nan

hp.write_map("sky_projection.fits", sky_map)
"""

ell = np.arange(0,LMAX + 500)
FWHM = np.deg2rad(8.0/60.0)
SIGMA = FWHM/2.35482
beam = np.exp(-0.5*ell*(ell + 1)*SIGMA**2)

hp.mwrfits("beam_file.fits", [beam])

bins = np.arange(0, LMAX, 20)
bins = bins.astype(float)
hp.mwrfits("bin_file.fits", [bins])

"""
output_folder = "/global/scratch2/sd/banerji/leap_output/2015-08-10--13-38-01_iqu_maps/"

weight_map = hp.read_map(output_folder + "weight_map.fits", field = 1)
hp.write_map("pol_weight_map.fits", weight_map)

hit_map = hp.read_map("nhits_map.fits")
binary_map = hit_map.copy()
binary_map[hit_map != 0] = 1
hp.write_map("binary_mask_map.fits", binary_map)
"""
