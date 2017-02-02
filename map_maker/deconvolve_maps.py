import numpy as np
import healpy as hp
import simulation.power_spectrum.spectra_tools as st
import os

sim_dir = "/global/homes/b/banerji/simulation/output/beam_set_CENTRAL"

map_dir_list = ["rec_1a1b", "rec_1a1b_sym"]

map_fwhm_list = np.radians([7.12, 7.12])/60.0

for i in range(len(map_dir_list)):
    print "Doing :", map_dir_list[i]
    sky_map_in = hp.read_map(os.path.join(sim_dir, map_dir_list[i], "sky_map.fits"))
    sky_map_in[np.isnan(sky_map_in)] = 0.0
    sky_map_out = st.deconvolve_map(sky_map_in, fwhm_in=map_fwhm_list[i], fwhm_out=0.0, lmax=3000)
    #deconvolve_map(map_in, fwhm_in=, fwhm_out=0.0, lmax=None, binary_mask=None, pol=False, wiener=True, sky_prior=None)
    hp.write_map(os.path.join(sim_dir, map_dir_list[i], "sky_map_deconvolved.fits"), sky_map_out)
