#!/usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from simulation.settings.global_settings import global_paths
import os

nside1 = 1024
nside2 = 4096
lmax = 4000

spectra = np.load(os.path.join(global_paths.spectra_folder, "test_4000.npy"))[0]
np.random.seed(1234)
map_low = hp.synfast(spectra, nside1, lmax=4000, pol=False, new=True)
np.random.seed(1234)
map_high = hp.synfast(spectra, nside2, lmax=4000, pol=False, new=True)

theta, phi = hp.pix2ang(nside2, np.arange(hp.nside2npix(nside2)))

map_comp = hp.get_interp_val(map_low, theta, phi)
