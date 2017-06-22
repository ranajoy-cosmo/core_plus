#!/usr/bin/env python

import numpy as np
import healpy as hp
import os
import sys
import simulation.power_spectrum.spectra_tools as st
from distutils.util import strtobool

sim_dir_name = sys.argv[1]
action = sys.argv[2]
if action == "create":
    beam_fwhm = float(sys.argv[3])
else:
    pass

scratch_dir = "/scratch1/scratchdirs/banerji" 
sim_dir = os.path.join(scratch_dir, "core_output", sim_dir_name)
mask_dir = os.path.join(scratch_dir, "core_maps")

verbose = False

sky_1a = hp.read_map(os.path.join(sim_dir, "rec_1a", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_1b = hp.read_map(os.path.join(sim_dir, "rec_1b", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_2a = hp.read_map(os.path.join(sim_dir, "rec_2a", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_2b = hp.read_map(os.path.join(sim_dir, "rec_2b", "sky_map.fits"), field=(0,1,2), verbose=verbose)
#sky_3a = hp.read_map(os.path.join(sim_dir, "rec_3a", "sky_map.fits"), field=(0,1,2), verbose=verbose)
#sky_3b = hp.read_map(os.path.join(sim_dir, "rec_3b", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_pair_1 = hp.read_map(os.path.join(sim_dir, "rec_pair_1", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_pair_2 = hp.read_map(os.path.join(sim_dir, "rec_pair_2", "sky_map.fits"), field=(0,1,2), verbose=verbose)
#sky_pair_3 = hp.read_map(os.path.join(sim_dir, "rec_pair_3", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_all = hp.read_map(os.path.join(sim_dir, "rec_all", "sky_map.fits"), field=(0,1,2), verbose=verbose)
sky_mask = hp.read_map(os.path.join(mask_dir, "mask25_1024.fits"), verbose=verbose) 

sky_avg_1 = np.empty((3, 12*1024**2))
sky_avg_2 = np.empty((3, 12*1024**2))
#sky_avg_3 = np.empty((3, 12*1024**2))
sky_avg_all = np.empty((3, 12*1024**2))
sky_leak_1 = np.empty((3, 12*1024**2))
sky_leak_2 = np.empty((3, 12*1024**2))
#sky_leak_3 = np.empty((3, 12*1024**2))
sky_leak_all = np.empty((3, 12*1024**2))

for i in range(3):
    sky_avg_1[i] = 0.5*(sky_1a[i] + sky_1b[i])
    sky_avg_2[i] = 0.5*(sky_2a[i] + sky_2b[i])
    #sky_avg_3[i] = 0.5*(sky_3a[i] + sky_3b[i])
    #sky_avg_all[i] = (1.0/6.0)*(sky_1a[i] + sky_1b[i] + sky_2a[i] + sky_2b[i] + sky_3a[i] + sky_3b[i])
    sky_avg_all[i] = 0.25*(sky_1a[i] + sky_1b[i] + sky_2a[i] + sky_2b[i])
    sky_leak_1[i] = sky_pair_1[i] - sky_avg_1[i]
    sky_leak_2[i] = sky_pair_2[i] - sky_avg_2[i]
    #sky_leak_3[i] = sky_pair_3[i] - sky_avg_3[i]
    sky_leak_all[i] = sky_all[i] - sky_avg_all[i]

del sky_1a
del sky_1b
del sky_2a
del sky_2b
#del sky_3a
#del sky_3b

if action=="load":
    spectra_leak_1 = np.load(os.path.join(sim_dir, "spectra_leak_1.npy"))
    spectra_leak_2 = np.load(os.path.join(sim_dir, "spectra_leak_2.npy"))
    #spectra_leak_3 = np.load(os.path.join(sim_dir, "spectra_leak_3.npy"))
    spectra_leak_all = np.load(os.path.join(sim_dir, "spectra_leak_all.npy"))
elif action=="create":
    spectra_leak_1 = st.estimate_cl(sky_map=sky_leak_1, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    spectra_leak_2 = st.estimate_cl(sky_map=sky_leak_2, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    #spectra_leak_3 = st.estimate_cl(sky_map=sky_leak_3, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    spectra_leak_all = st.estimate_cl(sky_map=sky_leak_all, lmax=2000, binary_mask=sky_mask, fwhm=np.radians(beam_fwhm/60.0), pol=True) 
    np.save(os.path.join(sim_dir, "spectra_leak_1"), spectra_leak_1)
    np.save(os.path.join(sim_dir, "spectra_leak_2"), spectra_leak_2)
    #np.save(os.path.join(sim_dir, "spectra_leak_3"), spectra_leak_3)
    np.save(os.path.join(sim_dir, "spectra_leak_all"), spectra_leak_all)
else:
    print "Provide correct action : (load/create)"
