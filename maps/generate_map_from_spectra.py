#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os
from simulation.lib.geometry.conversions import am2rad

def make_maps():

    # The input spectrum was generated using recent Planck parameters and using pycamb, the python wrapper around CAMB
    # The spectrum was generated from l=2 to l=4000
    # Thus we are assuming that the maps is created by healpy starting from l=2. ???Need to verify???
    input_spectra = np.load(map_params.spectra_file)
    map_seed = 1234
    
    if map_params.do_unbeamed_map:
        np.random.seed(map_seed)
        sky_map = hp.synfast(input_spectra, nside=map_params.nside, lmax=map_params.lmax, new=True, pol=map_params.do_pol, pixwin=map_params.do_pixel_window)
        sky_map = np.array(sky_map, dtype=np.float32)
        if map_params.write_map:
            map_file_name = map_params.map_file + str(map_params.nside) + '_0'
            if map_params.do_pixel_window:
                map_file_name = map_file_name + '_pixwin'
            hp.write_map(map_file_name + '.fits', sky_map)


    if map_params.do_beamed_map:
        np.random.seed(map_seed)
        sky_map_beamed = hp.synfast(input_spectra, nside=map_params.nside, lmax=map_params.lmax, fwhm=am2rad(map_params.fwhm), new=True, pol=map_params.do_pol, pixwin=map_params.do_pixel_window)
        sky_map_beamed = np.array(sky_map_beamed, dtype=np.float32)
        if map_params.write_map:
            map_file_name = map_params.map_file + str(map_params.nside) + '_' + str(int(map_params.fwhm))
            if map_params.do_pixel_window:
                map_file_name = map_file_name + '_pixwin'
            hp.write_map(map_file_name + '.fits', sky_map_beamed)


    if map_params.do_degraded_map:
        if map_params.do_unbeamed_map:
            sky_map_degraded = hp.ud_grade(sky_map, map_params.nside_deg)
            if map_params.write_map:
                map_file_name = map_params.map_file + str(map_params.nside_deg) + '_0_deg.fits'
                hp.write_map(map_file_name, sky_map_degraded)
        if map_params.do_beamed_map:
            sky_map_beamed_degraded = hp.ud_grade(sky_map_beamed, map_params.nside_deg)
            if map_params.write_map:
                map_file_name = map_params.map_file + str(map_params.nside_deg) + '_' + str(int(map_params.fwhm)) +"_deg.fits"
                hp.write_map(map_file_name, sky_map_beamed_degraded)
       
    
if __name__=='__main__':
    from custom_params import map_params
    make_maps()
