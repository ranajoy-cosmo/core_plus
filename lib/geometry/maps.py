#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

def change_system(map_in, coord_in, coord_out, pol=True):
    nside = hp.get_nside(map_in)
    out_pix = np.arange(12*nside**2)
    lat_out, lon_out = hp.pix2ang(nside, out_pix)
    lat_in, lon_in = rotate_lat_lon(lat_out, lon_out, coord_out, coord_in)
    if pol:
        map_out = np.empty((3, 12*nside**2))
        for i in range(3):
            map_out[i] = hp.get_interp_val(map_in[i], lat_in, lon_in)
        u_s_out, u_e_out = get_local_axes(lat_out, lon_out)
    else:
        map_out = hp.get_interp_val(map_in, lat_in, lon_in)
    return map_out

def rotate_vector(v_in, coord_in, coord_out):
    r = hp.Rotator(coord=[coord_in, coord_out])
    v_out = r(

def rotate_lat_lon(lat_in, lon_in, coord_in, coord_out):
    r = hp.Rotator(coord=[coord_in, coord_out])
    lat_out, lon_out = r(lat_in, lon_in)
    return lat_out, lon_out

def get_local_axes(lat, lon, degrees=False):
    if degrees:
        lat = np.radians(lat)
        lon = np.radians(lon)
    u_s = np.array(zip(np.cos(lat)*np.cos(lon), -np.cos(lat)*np.sin(lon), -np.sin(lat))) 
    u_e = np.array(zip(np.sin(lon), np.cos(lon), np.zeros(lon.size)))

    return u_s, u_e
