#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import simulation.power_spectrum.spectra_tools as st
from distutils.util import strtobool

sim_dir_name = sys.argv[1]
action = sys.argv[2]
noise = strtobool(sys.argv[3])

scratch_dir = "/scratch1/scratchdirs/banerji" 
sim_dir = os.path.join(scratch_dir, "core_output", sim_dir_name)
mask_dir = os.path.join(scratch_dir, "core_maps")

verbose = False

sky_a = hp.read_map(os.path.join(sim_dir, "rec_a", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_b = hp.read_map(os.path.join(sim_dir, "rec_b", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_template = hp.read_map(os.path.join(sim_dir, "rec_TEMPLATE", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_template_QU = hp.read_map(os.path.join(sim_dir, "rec_TEMPLATE_QU", "sky_map.fits"), field=(0,1), verbose=verbose)
sky_pair = hp.read_map(os.path.join(sim_dir, "rec_pair", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_diff = hp.read_map(os.path.join(sim_dir, "rec_diff_QU", "sky_map.fits"), field=(0,1), verbose=verbose)
sky_corr = hp.read_map(os.path.join(sim_dir, "rec_corrected", "sky_map.fits"), field=(0,1), verbose=verbose)
sky_mask = hp.read_map(os.path.join(mask_dir, "mask25_1024.fits"), verbose=verbose) 
if noise:
    sky_a_noise = hp.read_map(os.path.join(sim_dir, "rec_a_noise", "sky_map.fits"), field=(0,1,2), verbose=verbose)
    sky_b_noise = hp.read_map(os.path.join(sim_dir, "rec_b_noise", "sky_map.fits"), field=(0,1,2), verbose=verbose)
    sky_diff_noise = hp.read_map(os.path.join(sim_dir, "rec_diff_QU_noise", "sky_map.fits"), field=(0,1), verbose=verbose)
    sky_pair_noise = hp.read_map(os.path.join(sim_dir, "rec_pair_noise", "sky_map.fits"), field=(0,1,2), verbose=verbose)
est_y = np.load(os.path.join(sim_dir, "estimated_y.npy"))

sky_avg = np.empty((3, 12*1024**2))
sky_leak = np.empty((3, 12*1024**2))
sky_res_leak = np.zeros((3, 12*1024**2))

for i in range(3):
    sky_avg[i] = 0.5*(sky_a[i] + sky_b[i])
    if noise:
        sky_avg[i] -= 0.5*(sky_a_noise[i] + sky_b_noise[i])
    sky_leak[i] = sky_pair[i] - sky_avg[i]
    if noise:
        sky_leak[i] -= sky_pair_noise[i]

for i in range(2):
    sky_res_leak[i+1] = sky_corr[i] - sky_avg[i+1]
    if noise:
        sky_res_leak[i+1] -= sky_diff_noise[i]

if action=="load":
    spectra_leak = np.load(os.path.join(sim_dir, "spectra_leak.npy"))
    spectra_res_leak = np.load(os.path.join(sim_dir, "spectra_residual_leak.npy"))
elif action=="create":
    spectra_leak = st.estimate_cl(sky_map=sky_leak, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(7.68/60.0), pol=True) 
    spectra_res_leak = st.estimate_cl(sky_map=sky_res_leak, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(7.68/60.0), pol=True) 
    np.save(os.path.join(sim_dir, "spectra_leak"), spectra_leak)
    np.save(os.path.join(sim_dir, "spectra_residual_leak"), spectra_res_leak)
else:
    print "Provide correct action : (load/create)"
