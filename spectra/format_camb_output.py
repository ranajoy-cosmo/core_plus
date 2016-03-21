#!/usr/bin/env python 

import numpy as np

input_file = "r_001/test_totCls.dat"
out_file = "r_001/unlensed_cls"

read_spectra = np.loadtxt(input_file).T

ell = read_spectra[0]

norm_factor = ell*(ell+1)/(2*np.pi)
spectra = np.zeros((read_spectra.shape[0]-1, read_spectra.shape[1]+2))
spectra[...,2:] = read_spectra[1:]/norm_factor


np.save(out_file, spectra)
