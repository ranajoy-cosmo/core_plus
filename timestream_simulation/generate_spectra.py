import numpy as np
import healpy as hp
import simulation.power_spectrum.spectra_tools as st
import simulation.beam.fit_fwhm_to_spectra as ft
import sys
import os

sim_folder = "/global/homes/b/banerji/simulation/output/beam_set_CENTRAL"

sky_symm_folder = sys.argv[1]
sky_asym_folder = sys.argv[2]
sky_estm_folder = sys.argv[3]

sky_symm = hp.read_map(os.path.join(sim_folder, sky_symm_folder, "sky_map.fits"), field=(0,1,2))
sky_asym = hp.read_map(os.path.join(sim_folder, sky_asym_folder, "sky_map.fits"), field=(0,1,2))
sky_estm = hp.read_map(os.path.join(sim_folder, sky_estm_folder, "sky_map.fits"), field=(0,1,2))

sky_leakage = np.empty((3, 12*4096**2))
sky_res_leakage = np.empty((3, 12*4096**2))

for i in range(3):
    sky_leakage[i] = sky_asym[i] - sky_symm[i]
    sky_res_leakage[i] = sky_asym[i] - sky_symm[i] - sky_estm[i]

spectra_asym = st.estimate_cl(sky_asym, lmax=3000)

out = ft.estimate_fwhm(spectra_asym[0, 2:1501], 1500)

fwhm = np.radians(out.params['fwhm'].value/60.0)

spectra_asym = st.deconvolve_spectra(spectra_asym, fwhm, lmax=3000)
np.save(os.path.join(sim_folder, sky_asym_folder, "spectra_sky"), spectra_asym)

spectra_leakage = st.estimate_cl(sky_leakage, lmax=3000, fwhm=fwhm)
np.save(os.path.join(sim_folder, sky_asym_folder, "spectra_leakage"), spectra_leakage)

spectra_res_leakage = st.estimate_cl(sky_res_leakage, lmax=3000, fwhm=fwhm)
np.save(os.path.join(sim_folder, sky_estm_folder, "spectra_res_leakage"), spectra_res_leakage)
