#! /usr/bin/env python
import numpy as np
from simulation.lib.utilities.generic_class import Generic

settings = Generic()

settings.lmax = 2500

settings.nside = 1024

fwhm_arcmin = 8.0

settings.fwhm = np.deg2rad(fwhm_arcmin/60.0)

settings.sigma = settings.fwhm/2.35482

root_folder = "/Users/banerji/CORE+/simulation/"

input_spectra_folder = root_folder + "maps_and_spectra/spectra/"
settings.input_spectra_file = input_spectra_folder + "spectra_pure_E.npy"

output_map_folder = root_folder + "maps_and_spectra/maps/"
settings.output_map_file = output_map_folder + "test_map.fits"

settings.do_pol = True

settings.do_beam_smoothing = True
