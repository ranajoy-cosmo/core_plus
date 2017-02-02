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
if action=="create":
    beam_fwhm = float(sys.argv[4])

scratch_dir = "/scratch1/scratchdirs/banerji" 
sim_dir = os.path.join(scratch_dir, "core_output", sim_dir_name)
mask_dir = os.path.join(scratch_dir, "core_maps")

verbose = False

sky_1a = hp.read_map(os.path.join(sim_dir, "rec_1a", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_1b = hp.read_map(os.path.join(sim_dir, "rec_1b", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_2a = hp.read_map(os.path.join(sim_dir, "rec_2a", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_2b = hp.read_map(os.path.join(sim_dir, "rec_2b", "sky_map.fits"), field=(0,1,2), verbose=verbose)
#sky_template = hp.read_map(os.path.join(sim_dir, "rec_TEMPLATE_1", "sky_map.fits"), field=(0,1,2), verbose=verbose)
#sky_template_QU_1 = hp.read_map(os.path.join(sim_dir, "rec_TEMPLATE_QU_1", "sky_map.fits"), field=(0,1), verbose=verbose)
#sky_template_QU_2 = hp.read_map(os.path.join(sim_dir, "rec_TEMPLATE_QU_2", "sky_map.fits"), field=(0,1), verbose=verbose)
sky_pair_1 = hp.read_map(os.path.join(sim_dir, "rec_pair_1", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_pair_2 = hp.read_map(os.path.join(sim_dir, "rec_pair_2", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_quad = hp.read_map(os.path.join(sim_dir, "rec_quad", "sky_map.fits"), field=(0,1,2), verbose=verbose)
#sky_diff = hp.read_map(os.path.join(sim_dir, "rec_diff_QU", "sky_map.fits"), field=(0,1), verbose=verbose)
sky_corr_1 = hp.read_map(os.path.join(sim_dir, "rec_corrected_0001", "sky_map.fits"), field=(0,1), verbose=verbose)
sky_corr_2 = hp.read_map(os.path.join(sim_dir, "rec_corrected_0002", "sky_map.fits"), field=(0,1), verbose=verbose)
sky_corr = hp.read_map(os.path.join(sim_dir, "rec_corrected", "sky_map.fits"), field=(0,1), verbose=verbose)
sky_mask = hp.read_map(os.path.join(mask_dir, "mask25_1024.fits"), verbose=verbose) 

if noise:
    sky_a_noise = hp.read_map(os.path.join(sim_dir, "rec_a_noise", "sky_map.fits"), field=(0,1,2), verbose=verbose)
    sky_b_noise = hp.read_map(os.path.join(sim_dir, "rec_b_noise", "sky_map.fits"), field=(0,1,2), verbose=verbose)
    sky_diff_noise = hp.read_map(os.path.join(sim_dir, "rec_diff_QU_noise", "sky_map.fits"), field=(0,1), verbose=verbose)
    sky_pair_noise = hp.read_map(os.path.join(sim_dir, "rec_pair_noise", "sky_map.fits"), field=(0,1,2), verbose=verbose)
#est_y = np.load(os.path.join(sim_dir, "estimated_y.npy"))

sky_avg_1 = np.empty((3, 12*1024**2))
sky_avg_2 = np.empty((3, 12*1024**2))
sky_avg = np.empty((3, 12*1024**2))
sky_leak_1 = np.empty((3, 12*1024**2))
sky_leak_2 = np.empty((3, 12*1024**2))
sky_leak = np.empty((3, 12*1024**2))
sky_res_leak_1 = np.zeros((3, 12*1024**2))
sky_res_leak_2 = np.zeros((3, 12*1024**2))
sky_res_leak = np.zeros((3, 12*1024**2))

for i in range(3):
    sky_avg_1[i] = 0.5*(sky_1a[i] + sky_1b[i])
    sky_avg_2[i] = 0.5*(sky_2a[i] + sky_2b[i])
    sky_avg[i] = 0.5*(sky_avg_2[i] + sky_avg_1[i])
    if noise:
        sky_avg[i] -= 0.5*(sky_a_noise[i] + sky_b_noise[i])
    sky_leak_1[i] = sky_pair_1[i] - sky_avg_1[i]
    sky_leak_2[i] = sky_pair_2[i] - sky_avg_2[i]
    sky_leak[i] = sky_quad[i] - sky_avg[i]
    if noise:
        sky_leak[i] -= sky_pair_noise[i]

for i in range(2):
    sky_res_leak_1[i+1] = sky_corr_1[i] - sky_avg_1[i+1]
    sky_res_leak_2[i+1] = sky_corr_2[i] - sky_avg_2[i+1]
    sky_res_leak[i+1] = sky_corr[i] - sky_avg[i+1]
    if noise:
        sky_res_leak[i+1] -= sky_diff_noise[i]

if action=="load":
    spectra_leak = np.load(os.path.join(sim_dir, "spectra_leak.npy"))
    spectra_res_leak = np.load(os.path.join(sim_dir, "spectra_residual_leak.npy"))
elif action=="create":
    spectra_leak_1 = st.estimate_cl(sky_map=sky_leak_1, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    spectra_leak_2 = st.estimate_cl(sky_map=sky_leak_2, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    spectra_leak = st.estimate_cl(sky_map=sky_leak, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    spectra_res_leak_1 = st.estimate_cl(sky_map=sky_res_leak_1, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    spectra_res_leak_2 = st.estimate_cl(sky_map=sky_res_leak_2, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    spectra_res_leak = st.estimate_cl(sky_map=sky_res_leak, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    np.save(os.path.join(sim_dir, "spectra_leak_pair1"), spectra_leak)
    np.save(os.path.join(sim_dir, "spectra_leak_pair2"), spectra_leak)
    np.save(os.path.join(sim_dir, "spectra_leak"), spectra_leak)
    np.save(os.path.join(sim_dir, "spectra_residual_leak_pair1"), spectra_res_leak)
    np.save(os.path.join(sim_dir, "spectra_residual_leak_pair2"), spectra_res_leak)
    np.save(os.path.join(sim_dir, "spectra_residual_leak"), spectra_res_leak)
    pass
else:
    print "Provide correct action : (load/create)"
