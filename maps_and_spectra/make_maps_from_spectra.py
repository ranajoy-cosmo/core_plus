#! /usr/bin/env python 

import numpy as np
import healpy as hp

from map_settings import settings

input_spectra = np.load(settings.input_spectra_file)

if settings.do_beam_smoothing:
    sky_map = hp.synfast(input_spectra, nside = settings.nside, lmax = settings.lmax, fwhm = settings.fwhm, new = True, 
                pol = settings.do_pol)
else:
    sky_map = hp.synfast(input_spectra, nside = settings.nside, lmax = settings.lmax, new = True, pol = settings.do_pol)

hp.write_map(settings.output_map_file, sky_map)


