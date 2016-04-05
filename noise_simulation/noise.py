#!/usr/bin/env python 

import numpy as np
import healpy as hp

def noise_to_map(nside, noise_level, pol):
    pix_size = hp.nside2resol(nside, arcmin=True)
    noise_sigma = noise_level/np.sqrt(pix_size)
    npix = hp.nside2npix(nside)
    if pol:
        noise_map = np.random.normal(loc=0.0, scale=noise_sigma, size=(3, npix))
    else:
        noise_map = np.random.normal(loc=0.0, scale=noise_sigma, size=npix)
    return noise_map

def noise_to_timestream(length, white_noise_sigma, f_knee=None, alpha=None, type="white"):
    if type is "white":
        noise_ts = np.random.normal(loc=0.0, scale=white_noise_sigma, size=length)

    return noise_ts
