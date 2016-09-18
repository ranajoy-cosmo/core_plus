import numpy as np
import healpy as hp
import os
from simulation.params.custom_params import global_paths

OMEGA_SKY = 4*np.pi                         #radian^2
OMEGA_DEG_SQ = np.radians(1.0)**2           #radian^2
OMEGA_ARCMIN_SQ = np.radians(1.0/60.0)**2   #radian^2

factor = 2*np.sqrt(2*np.log(2))                                     #FWHM <--> sigma

spectra_folder = os.path.join(global_paths.base_dir, "spectra")
spectra_file = os.path.join(spectra_folder, "r_001/lensedtot_cls.npy")

def get_pixel_solid_angle(nside, arcmin=False):
    pix_res = hp.nside2res(nside, arcmin=True)
    if arcmin:
        return pix_res**2
    else:
        return OMEGA_ARCMIN_SQ*pix_res**2

def time_per_pixel(nside, t_mission):
    omega_pix = get_pixel_solid_angle(nside, arcmin=True)
    t_pix = t_mission*omega_pix/OMEGA_SKY
    return t_pix

def time_per_arcmin_sq(t_mission):
    t_arcmin_sq = t_mission*OMEGA_ARCMIN_SQ/OMEGA_SKY
    return t_arcmin_sq

def sensitivity_to_noise_per_pix(det_sens, t_mission, nside, n_det=1):
    t_pix = time_per_pixel(nside, t_mission)
    noise_per_pix = det_sens / np.sqrt(t_pix) / np.sqrt(n_det)
    return noise_per_pix

def sensitivity_to_noise_arcmin_sq(det_sens, t_mission, n_det=1):
    t_arcmin_sq = time_per_arcmin_sq(t_mission)
    noise_arcmin_sq = det_sen / np.sqrt(t_arcmin_sq) / np.sqrt(n_det)
    return noise_arcmin_sq

def get_noise_per_solid_angle(det_sens, t_mission, n_det=1):
    noise_per_solid_angle = det_sens * np.sqrt(OMEGA_SKY/t_mission) / np.sqrt(n_det)
    return noise_per_solid_angle

def get_sensitivity_from_noise_arcmin_sq(noise_arcmin_sq, t_mission, n_det=1):
    t_arcmin_sq = time_per_arcmin_sq(t_mission)
    det_sens = noise_arcmin_sq * np.sqrt(t_arcmin_sq) * np.sqrt(n_det)
    return det_sens

def get_variance_Cl(det_sens, beam_fwhm, lmax):
    Cl = np.load(spectra_file)[..., :lmax+1]
    ell = np.arange(lmax+1)
    w = 1.0 / get_noise_per_solid_angle**2
    Bl = hp.gauss_beam(fwhm=beam_fwhm, lmax=lmax, pol=True)
    
    variance_TT = (2/(2*ell+1))*(Cl[0] + Bl**2/w)
    return variance_TT

def get_spectral_noise_per_multipole(set_sens, beam_fwhm, lmax):
    Cl = np.load(spectra_file)[..., :lmax+1]
    variance_TT = get_variance_Cl(det_sens, beam_fwhm , lmax)
    spectra_noise = np.sqrt(varaince_TT) / Cl[0] 

    return spectral_noise
