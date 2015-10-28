#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os

def make_maps(settings=None):
    if settings is None:
        from custom_settings import settings
    input_spectra = np.load(settings.spectra_file)
    map_seed = 1234
    
    if settings.do_unbeamed_map:
        np.random.seed(map_seed)
        sky_map = hp.synfast(input_spectra, nside=settings.nside, lmax=3950, new=True, pol=settings.do_pol)
        sky_map = np.array(sky_map, dtype=np.float32)

    if settings.do_beamed_map:
        np.random.seed(map_seed)
        sky_map_beamed = hp.synfast(input_spectra, nside=settings.nside, lmax=3950, fwhm=settings.fwhm, new=True, pol=settings.do_pol)
        sky_map_beamed = np.array(sky_map_beamed, dtype=np.float32)

    if settings.do_degraded_map:
        sky_map_degraded = hp.ud_grade(sky_map, nside_deg)
        sky_map_beamed_degraded = hp.ud_grade(sky_map_beamed, nside_deg)
       
    if settings.write_map:
        if settings.do_unbeamed_map:
            hp.write_map(settings.map_file + '_0.fits', sky_map)
        if settings.do_beamed_map:
            hp.write_map(settings.map_file + '_' + str(int(settings.fwhm_arcmin)) +".fits", sky_map_beamed)
        if settings.do_unbeamed_map and settings.do_degraded_map:
            hp.write_map(settings.map_file + '_deg_0.fits', sky_map_degraded)
        if settings.do_beamed_map and settings.do_degraded_map:
            hp.write_map(settings.map_file + '_deg_' + str(int(settings.fwhm_arcmin)) +".fits", sky_map_beamed_degraded)
    
    if settings.display_maps:
        if settings.do_pol:
            for i in (0,1,2):
                hp.mollzoom(sky_map[i])
                plt.show()
        else:
            hp.mollzoom(sky_map)
            plt.show()
    
    if settings.return_map:
        return sky_map

if __name__=='__main__':
    from custom_settings import settings
    make_maps(settings)
