#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os

def make_maps(settings=None):
    if settings is None:
        from custom_settings import settings

    # The input spectrum was generated using recent Planck parameters and using pycamb, the python wrapper around CAMB
    # The spectrum was generated from l=2 to l=4000
    # Thus we are assuming that the maps is created by healpy starting from l=2. ???Need to verify???
    input_spectra = np.load(settings.spectra_file)
    map_seed = 1234
    
    if settings.do_unbeamed_map:
        np.random.seed(map_seed)
        sky_map = hp.synfast(input_spectra, nside=settings.nside, lmax=settings.lmax, new=True, pol=settings.do_pol)
        sky_map = np.array(sky_map, dtype=np.float32)
        if settings.write_map:
            map_file_name = settings.map_file + str(settings.nside) + '_0.fits'
            hp.write_map(map_file_name, sky_map)


    if settings.do_beamed_map:
        np.random.seed(map_seed)
        sky_map_beamed = hp.synfast(input_spectra, nside=settings.nside, lmax=settings.lmax, fwhm=settings.fwhm, new=True, pol=settings.do_pol)
        sky_map_beamed = np.array(sky_map_beamed, dtype=np.float32)
        if settings.write_map:
            map_file_name = settings.map_file + str(settings.nside) + '_' + str(int(settings.fwhm_arcmin)) +".fits"
            hp.write_map(map_file_name, sky_map_beamed)


    if settings.do_degraded_map:
        if settings.do_unbeamed_map:
            sky_map_degraded = hp.ud_grade(sky_map, settings.nside_deg)
            if settings.write_map:
                map_file_name = settings.map_file + str(settings.nside_deg) + '_0_deg.fits'
                hp.write_map(map_file_name, sky_map_degraded)
        if settings.do_beamed_map:
            sky_map_beamed_degraded = hp.ud_grade(sky_map_beamed, settings.nside_deg)
            if settings.write_map:
                map_file_name = settings.map_file + str(settings.nside_deg) + '_' + str(int(settings.fwhm_arcmin)) +"_deg.fits"
                hp.write_map(map_file_name, sky_map_beamed_degraded)
       
    
if __name__=='__main__':
    from custom_settings import settings
    make_maps(settings)
