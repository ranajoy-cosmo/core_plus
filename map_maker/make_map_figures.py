#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import simulation.power_spectrum.spectra_tools as st
import simulation.power_spectrum.noise_power_spectra as nps
import simulation.lib.plotting.plot_th as pt
import sys
import os

sim_dir = "/global/homes/b/banerji/simulation/output/focal_plane_noise" 
map_folder_list = []
det_number_list = []
scan_time_list = []

general_plot_dir = "/global/homes/b/banerji/simulation/mission_paper_plots"

sky_map_title = ['I', 'Q', 'U']
cov_map_title = ['II', 'IQ', 'IU', 'QQ', 'QU', 'UU']

current_dir = "/global/homes/b/banerji/simulation/map_maker"

det_sens = 50.0
scan_rate = 85.0
nside = 1024
year = 365.25*24*60*60.0
precession = 4*24*60*60.0

noise_rms = noise_sens*np.sqrt(scan_rate)
pix_area = hp.nside2resol(nside, arcmin=True)**2

for map_folder in map_folder_list:

    map_dir = os.path.join(sim_dir, map_folder)
    plot_dir = os.path.join(general_plot_dir, map_folder)
    
    if os.path.exists(plot_dir):
        shutil.rmtree(plot_dir)
    os.makedirs(plot_dir)

    noise_map = hp.read_map(os.path.join(map_dir, "sky_map.fits"), field=(0,1,2))
    cov_map = hp.read_map(os.path.join(map_dir, "covariance_maps.fits"), field=(0,1,2,3,4,5))
    hitmap = hp.read_map(os.path.join(map_dir, "hitmap.fits"))

    mask = hitmap > 3

    spectra_obs = st.estimate_cl(noise_map, lmax=2000, binary_mask=hitmap>3, fwhm=np.radians(7.7/60.0))
    spectra_th = nps.sensitivity_to_noise_mean_Cl(det_sens, np.radians(7.7/60.0), 2000, t_mission=scan_time_dict[map_folder], n_det=det_number_dict[map_folder])
