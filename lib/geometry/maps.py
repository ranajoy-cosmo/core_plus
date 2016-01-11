#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

def rotate_map(map_in, angles, convention="ZYX"):
    nside = hp.get_nside(map_in)
    lat_in, lon_in = hp.pix2ang(nside, np.arange(12*nside**2))
    r = hp.Rotator(rot=angles, deg=True)
    lat_out, lon_out = r(lat_in, lon_in)
    pix_new = hp.ang2pix(nside, lat_out, lon_out)
    map_out = np.empty(map_in.size)
    map_out[pix_new] = map_in
    return map_out
