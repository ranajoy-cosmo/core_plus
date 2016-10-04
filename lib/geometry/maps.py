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

def change_system(map_in, coord_in, coord_out, pol=True):
    nside = hp.get_nside(map_in)
    out_pix = np.arange(12*nside**2)
    lat_out, lon_out = hp.pix2ang(nside, out_pix)
    r = hp.Rotator(coord=[coord_out, coord_in])
    lat_in, lon_in = r(lat_out, lon_out)
    if pol:
        map_out = np.empty((3, 12*nside**2))
        for i in range(3):
            map_out[i] = hp.get_interp_val(map_in[i], lat_in, lon_in)
    else:
        map_out = hp.get_interp_val(map_in, lat_in, lon_in)
    return map_out
