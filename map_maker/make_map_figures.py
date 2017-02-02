#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
import os

sim_dir = sys.argv[1]
map_folder_list = os.listdir(sim_dir)
#map_folder_list = ['all_bolos_one_precession']

sky_map_title = ['I', 'Q', 'U']
cov_map_title = ['II', 'IQ', 'IU', 'QQ', 'QU', 'UU']

current_dir = "/global/homes/b/banerji/simulation/map_maker"

for map_folder in map_folder_list:

    general_noise_map_folder = os.path.join(current_dir, "fp_noise_maps")
    general_cov_map_folder = os.path.join(current_dir, "fp_cov_maps")
    noise_map_folder = os.path.join(general_noise_map_folder, map_folder)
    cov_map_folder = os.path.join(general_cov_map_folder, map_folder)


    if not os.path.exists(general_noise_map_folder):
        os.makedirs(general_noise_map_folder)
    if not os.path.exists(noise_map_folder):
        os.makedirs(noise_map_folder)
    if not os.path.exists(general_cov_map_folder):
        os.makedirs(general_cov_map_folder)
    if not os.path.exists(cov_map_folder):
        os.makedirs(cov_map_folder)

#    print map_folder
#    print noise_map_folder
#    print cov_map_folder
#    print os.path.join(sim_dir, map_folder, "sky_map.fits")
    hitmap = hp.read_map(os.path.join(sim_dir, map_folder, "hitmap.fits"))
    noise_map = hp.read_map(os.path.join(sim_dir, map_folder, "sky_map.fits"), field=(0,1,2))
    cov_matrix_maps = hp.read_map(os.path.join(sim_dir, map_folder, "covariance_maps.fits"), field=(0,1,2,3,4,5))

    for i in range(3):
        plt.figure()
        noise_map[i][hitmap==0] = np.nan
        hp.mollview(noise_map[i], title=sky_map_title[i], unit='$\mu K$')
        plt.savefig(os.path.join(noise_map_folder, sky_map_title[i]))

    for i in range(6):
        plt.figure()
        cov_matrix_maps[i][hitmap==0] = np.nan
        #cov_matrix_maps[i][cov_matrix_maps[i]==1] = np.nan
        hp.mollview(cov_matrix_maps[i], title=cov_map_title[i])
        plt.savefig(os.path.join(cov_map_folder, cov_map_title[i]))
