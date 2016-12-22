import numpy as np
from lmfit import minimize, Parameters, fit_report
import simulation.power_spectrum.noise_power_spectra as nps
from simulation.global_config import global_paths
import sys
import os

spectra_out_file = sys.argv[1]
fid_spectra_folder = os.path.join(global_paths.base_dir, "spectra")
fid_spectra_file = os.path.join(fid_spectra_folder, "r_0001/lensedtot_cls.npy")              #Fiducial spectra
lmax = 2000

def residual(params, ell, out_spectra, var_spectra, fid_spectra):
    factor = 2*np.sqrt(2*np.log(2))
    beam_sigma = np.radians(params['fwhm']/60.0)/factor
    model = fid_spectra*np.exp(-ell*(ell+1)*beam_sigma**2)

    return (out_spectra - model)/np.sqrt(var_spectra)

ell = np.arange(lmax+1)[2:]
spectra_var = nps.cosmic_variance_Cl(lmax)[0,2:]
spectra_fid = np.load(fid_spectra_file)[0,2:lmax+1]
spectra_out = np.load(spectra_out_file)[2:]

#print ell.shape, spectra_var.shape, spectra_fid.shape, spectra_out.shape
#sys.exit()

params = Parameters()
params.add('fwhm', value=5.0)

out = minimize(residual, params, args=(ell, spectra_out, spectra_var, spectra_fid))

print(fit_report(out))
