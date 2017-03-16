#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import simulation.power_spectrum.spectra_tools as st

sim_dir = "/global/homes/b/banerji/simulation/output/focal_plane_noise" 

#map_dir_list = ['2_bolo_half_year', '36_bolo_1_precession', '36_bolo_1_year', '9_bolo_half_year', '36_bolo_half_year', '36_bolo_4_year', '9_bolo_1_year', '18_bolo_1_precession', '18_bolo_half_year', '1_bolo_1_year', '18_bolo_4_year', '9_bolo_4_year',  '18_bolo_1_year', '1_bolo_half_year', '2_bolo_1_year', '2_bolo_4_year', '4_bolo_1_year', '1_bolo_4_year']
#map_dir_list = ['2_bolo_half_year', '36_bolo_1_precession', '36_bolo_1_year', '9_bolo_half_year', '36_bolo_half_year', '9_bolo_1_year', '18_bolo_1_precession', '18_bolo_half_year', '1_bolo_1_year', '18_bolo_4_year', '9_bolo_4_year',  '18_bolo_1_year', '1_bolo_half_year', '2_bolo_1_year', '2_bolo_4_year', '4_bolo_1_year', '1_bolo_4_year']
map_dir_list = ['36_bolo_4_year']

for map_dir in map_dir_list:
    print "Doing {}\n".format(map_dir)
    map_dir_dir = os.path.join(sim_dir, map_dir)
    noise_map = hp.read_map(os.path.join(map_dir_dir, "sky_map.fits"), field=(0,1,2))
    hitmap = hp.read_map(os.path.join(map_dir_dir, "hitmap.fits"))

    spectra_obv = st.estimate_cl(noise_map, lmax=3000, binary_mask=hitmap>3, fwhm=np.radians(7.7/60.0))
    spectra_obv[..., :2] = 0.0
    
    np.save(os.path.join(map_dir_dir, "spectra_noise"), spectra_obv)
