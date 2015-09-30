#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

def make_maps(settings=None):
    input_spectra = np.load(settings.spectra_file)

    if settings.do_beam_smoothing:
        sky_map = hp.synfast(input_spectra, nside = settings.nside, lmax = settings.lmax, fwhm = settings.fwhm, new = True, 
                    pol = settings.do_pol)
    else:
        sky_map = hp.synfast(input_spectra, nside = settings.nside, lmax = settings.lmax, new = True, pol = settings.do_pol)

    if settings.write_map:
        hp.write_map(settings.map_file, sky_map)

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
    from settings_map import settings
    make_maps(settings)
