import numpy as np
import healpy as hp
import os
from simulation.global_config import global_paths

OMEGA_SKY = 4*np.pi                         #radian^2
OMEGA_ARCMIN_SQ = np.radians(1.0/60.0)**2   #radian^2

FWHM_FACTOR = 2*np.sqrt(2*np.log(2))                                     #FWHM = FWHM_FACTOR*sigma

T_YEAR = 365*24*60*60.0                     #seconds
N_DET_145 = 144

spectra_folder = os.path.join(global_paths.base_dir, "spectra")
spectra_file = os.path.join(spectra_folder, "r_0001/unlensed_cls.npy")


def pixel_solid_angle(nside, arcmin=False):
    pix_res = hp.nside2resol(nside, arcmin=True)
    if arcmin:
        pix_solid_angle = pix_res**2
    else:
        pix_solid_angle = OMEGA_ARCMIN_SQ*pix_res**2
    return pix_solid_angle


def time_per_pixel(nside, t_mission=T_YEAR):
    omega_pix = pixel_solid_angle(nside)
    t_pix = t_mission*omega_pix/OMEGA_SKY
    return t_pix


def time_per_arcmin_sq(t_mission=T_YEAR):
    t_arcmin_sq = t_mission*OMEGA_ARCMIN_SQ/OMEGA_SKY
    return t_arcmin_sq


def time_per_solid_angle(t_mission=T_YEAR):
    t_unit_solid_angle = t_mission/OMEGA_SKY
    return t_unit_solid_angle


def sensitivity_to_noise_pix(det_sens, nside, t_mission=T_YEAR, n_det=N_DET_145):
    t_pix = time_per_pixel(nside, t_mission)
    noise_pix_T = det_sens / np.sqrt(t_pix) / np.sqrt(n_det)
    factor = np.array([1.0, np.sqrt(2), np.sqrt(2)])
    return noise_pix_T*factor


def sensitivity_to_noise_arcmin(det_sens, t_mission=T_YEAR, n_det=N_DET_145):
    t_arcmin_sq = time_per_arcmin_sq(t_mission)
    noise_arcmin_T = det_sens / np.sqrt(t_arcmin_sq) / np.sqrt(n_det)
    factor = np.array([1.0, np.sqrt(2), np.sqrt(2)])
    return noise_arcmin_T*factor


def sensitivity_to_noise_solid_angle(det_sens, t_mission=T_YEAR, n_det=N_DET_145):
    t_unit_solid_angle = time_per_solid_angle(t_mission)
    noise_solid_angle_T = det_sens / np.sqrt(t_unit_solid_angle) / np.sqrt(n_det)
    factor = np.array([1.0, np.sqrt(2), np.sqrt(2)])
    return noise_solid_angle_T*factor


def noise_arcmin_to_noise_pix(noise_arcmin, nside): 
    if not type(noise_arcmin) == np.ndarray:
        noise_arcmin = np.array([noise_arcmin])
    omega_pix = pixel_solid_angle(nside)
    noise_pix = noise_arcmin*np.sqrt(OMEGA_ARCMIN_SQ/omega_pix)
    return noise_pix


def noise_arcmin_to_noise_solid_angle(noise_arcmin):
    if not type(noise_arcmin) == np.ndarray:
        noise_arcmin = np.array([noise_arcmin])
    noise_solid_angle = noise_arcmin*np.sqrt(OMEGA_ARCMIN_SQ)
    return noise_solid_angle


def noise_arcmin_to_sensitivity(noise_arcmin, t_mission=T_YEAR, n_det=N_DET_145):
    t_arcmin_sq = time_per_arcmin_sq(t_mission)
    det_sens = noise_arcmin * np.sqrt(t_arcmin_sq) * np.sqrt(n_det)
    return det_sens


def sensitivity_to_variance_Cl(det_sens, beam_fwhm, lmax, t_mission=T_YEAR, n_det=N_DET_145, f_sky=1.0, spectra_fid=None):
    if spectra_fid == None:
        spectra_fid = spectra_file
    Cl = np.load(spectra_fid)[..., :lmax+1]
    ell = np.arange(lmax+1)
    w_inv = sensitivity_to_noise_solid_angle(det_sens, t_mission, n_det)**2
    Bl = hp.gauss_beam(fwhm=beam_fwhm, lmax=lmax, pol=True)
    variance = np.empty((3, lmax+1))
    for i in range(3):
        variance[i] = (2.0/(2.0*ell+1)/f_sky)*(Cl[i] + w_inv[i]/Bl[...,i]**2)**2
    return variance


def sensitivity_to_noise_power_spectra(det_sens, beam_fwhm, lmax, t_mission=T_YEAR, n_det=N_DET_145, f_sky=1.0):
    ell = np.arange(lmax+1)
    w_inv = sensitivity_to_noise_solid_angle(det_sens, t_mission, n_det)**2
    Bl = hp.gauss_beam(fwhm=beam_fwhm, lmax=lmax, pol=True)
    variance = np.empty((3, lmax+1))
    for i in range(3):
        variance[i] = (2.0/(2.0*ell+1)/f_sky)*(w_inv[i]/Bl[...,i]**2)**2
    return variance


def sensitivity_to_fractional_noise_per_multipole(det_sens, beam_fwhm, lmax, t_mission=T_YEAR, n_det=N_DET_145, f_sky=1.0, spectra_fid=None):
    if spectra_fid == None:
        spectra_fid = spectra_file
    Cl = np.load(spectra_fid)[..., :lmax+1]
    variance = sensitivity_to_variance_Cl(det_sens, beam_fwhm , lmax, t_mission, n_det, f_sky, spectra_fid)
    fractional_noise = np.empty(variance.shape)
    for i in range(3):
        fractional_noise[i] = np.sqrt(variance[i]) / Cl[i] 
    return fractional_noise


def noise_arcmin_to_variance_Cl(noise_arcmin, beam_fwhm, lmax, f_sky=1.0, spectra_fid=None):
    if spectra_fid == None:
        spectra_fid = spectra_file
    Cl = np.load(spectra_fid)[..., :lmax+1]
    ell = np.arange(lmax+1)
    w_inv = noise_arcmin_to_noise_solid_angle(noise_arcmin)**2
    Bl = hp.gauss_beam(fwhm=beam_fwhm, lmax=lmax, pol=True)
    variance = np.empty((3, lmax+1))
    for i in range(3):
        variance[i] = (2.0/(2.0*ell+1)/f_sky)*(Cl[i] + w_inv[i]/Bl[...,i]**2)**2
    return variance


def noise_arcmin_to_fractional_noise_per_multipole(noise_arcmin, beam_fwhm, lmax, f_sky=1.0, spectra_fid=None):
    if spectra_fid == None:
        spectra_fid = spectra_file
    Cl = np.load(spectra_fid)[..., :lmax+1]
    variance = noise_arcmin_to_variance_Cl(noise_arcmin, beam_fwhm, lmax, f_sky, spectra_fid)
    fractional_noise = np.empty(variance.shape)
    for i in range(3):
        fractional_noise[i] = np.sqrt(variance[i]) / Cl[i] 
    return fractional_noise


def noise_arcmin_to_noise_maps(noise_arcmin, nside, pol=True):
    npix = 12*nside**2
    noise_pix = noise_arcmin_to_noise_pix(noise_arcmin, nside)
    if pol:
        sky = np.empty((3, npix))
        for i in range(3):
            sky[i] = np.random.normal(loc=0.0, scale=noise_pix[i], size=npix)  
    else:
        sky = np.random.normal(loc=0.0, scale=noise_pix[0], size=npix)  
    return sky
