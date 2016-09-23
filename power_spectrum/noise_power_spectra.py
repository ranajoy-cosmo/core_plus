import numpy as np
import healpy as hp
import os
from simulation.params.custom_params import global_paths

OMEGA_SKY = 4*np.pi                         #radian^2
OMEGA_ARCMIN_SQ = np.radians(1.0/60.0)**2   #radian^2

factor = 2*np.sqrt(2*np.log(2))                                     #FWHM <--> sigma

spectra_folder = os.path.join(global_paths.base_dir, "spectra")
spectra_file = os.path.join(spectra_folder, "r_001/lensedtot_cls.npy")

def get_pixel_solid_angle(nside, arcmin=False):
    pix_res = hp.nside2resol(nside, arcmin=True)
    if arcmin:
        return pix_res**2
    else:
        return OMEGA_ARCMIN_SQ*pix_res**2

def time_per_pixel(nside, t_mission):
    omega_pix = get_pixel_solid_angle(nside, arcmin=False)
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
    noise_arcmin_sq = det_sens / np.sqrt(t_arcmin_sq) / np.sqrt(n_det)
    return noise_arcmin_sq

def get_noise_per_solid_angle(det_sens, t_mission, n_det=1):
    noise_per_solid_angle = det_sens * np.sqrt(OMEGA_SKY/t_mission) / np.sqrt(n_det)
    return noise_per_solid_angle

def get_sensitivity_from_noise_arcmin_sq(noise_arcmin_sq, t_mission, n_det=1):
    t_arcmin_sq = time_per_arcmin_sq(t_mission)
    det_sens = noise_arcmin_sq * np.sqrt(t_arcmin_sq) * np.sqrt(n_det)
    return det_sens

def get_variance_Cl(det_sens, beam_fwhm, lmax, t_mission=365*24*60*60, n_det=1):
    Cl = np.load(spectra_file)[..., :lmax+1]
    ell = np.arange(lmax+1)
    w = 1.0 / get_noise_per_solid_angle(det_sens, t_mission, n_det)**2
    Bl = hp.gauss_beam(fwhm=beam_fwhm, lmax=lmax, pol=True)
    
    variance_TT = (2.0/(2.0*ell+1))*(Cl[0] + Bl[...,0]**2/w)
    return variance_TT

def get_spectral_noise_per_multipole(det_sens, beam_fwhm, lmax, t_mission, n_det=1):
    Cl = np.load(spectra_file)[..., :lmax+1]
    variance_TT = get_variance_Cl(det_sens, beam_fwhm , lmax, t_mission, n_det)
    spectral_noise = np.sqrt(variance_TT) / Cl[0] 

    return spectral_noise

def get_noise_maps(map_noise_level, nside, pol=True):
    npix = 12*nside**2
    sigma_pix = map_noise_level / np.sqrt(get_pixel_solid_angle(nside, arcmin=True))
    if pol:
        sky = np.empty((3, npix))
        sky[0] = np.random.normal(loc=0.0, scale=sigma_pix, size=npix)  
        sky[1] = np.random.normal(loc=0.0, scale=sigma_pix/np.sqrt(2), size=npix)  
        sky[2] = np.random.normal(loc=0.0, scale=sigma_pix/np.sqrt(2), size=npix)  
    else:
        sky = np.random.normal(loc=0.0, scale=map_noise_level, size=npix)  

    return sky

def get_noise_maps_from_sensitivity(det_sens, nside, t_mission, n_det=1, pol=True):
    sigma_arcmin = sensitivity_to_noise_arcmin_sq(det_sens, t_mission, n_det)
    return get_noise_maps(sigma_arcmin, nside, pol)
