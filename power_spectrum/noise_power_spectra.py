import numpy as np
import healpy as hp
import os
from simulation.global_config import global_paths

"""
This module calculates and converts the noise parameters for a CMB experiment and also generates the noise power spectra given such parameters.
Calculations are done assuming the observer is at the centre of a 2-Sphere with the last scattering surface being the surface of the sphere and at radius unity away from the observer.
Features include :
    1) Calculate the solid angle of a pixel of a Healpix map from the nside
    2) Time spent by the detector per pixel of the Healpix map for the given nside and time span
    3) Time spent by the detector per arcminute square of the sky for the given time span
    4) Time spent by the detector per unit solid angle of the sky for the given time span
    5) Detector sensitivity to noise rms in each pixel of the Healpix map for the given nside, time span and number of detectors
    6) Detector sensitivity to noise rms in uK.arcmin for the given time span and number of detectors
    7) Detector sensitivity to noise rms in uK.steradian for the given time span and number of detectors
"""
OMEGA_SKY = 4*np.pi                         #steradian
OMEGA_ARCMIN_SQ = np.radians(1.0/60.0)**2   #steradian

FWHM_FACTOR = 2*np.sqrt(2*np.log(2))                                     #FWHM = FWHM_FACTOR*sigma

T_MISSION = 4*365*24*60*60.0                     #seconds : 4 years
N_DET_145 = 144                                     #for 145GHz

spectra_folder = os.path.join(global_paths.base_dir, "spectra")
spectra_file = os.path.join(spectra_folder, "r_0001/lensedtot_cls.npy")              #Fiducial spectra

def pixel_solid_angle(nside, arcmin=False):
    """
    Gives the solid angle subtended by a pixel on a Healpix map with a given nside.
    Default units ; steradian
    if arcmin is True, units = arcmin^2
    """
    pix_solid_angle = hp.nside2resol(nside, arcmin)**2
#    pix_res = hp.nside2resol(nside, arcmin=True)
#    if arcmin:
#        pix_solid_angle = pix_res**2
#    else:
#        pix_solid_angle = OMEGA_ARCMIN_SQ*pix_res**2
    return pix_solid_angle


def time_per_pixel(nside, t_mission=T_MISSION, n_det=1):
    """
    Time spent by a detector per pixel during a period t_mission and assuming an uniform scan.
    Default units : seconds
    """
    omega_pix = pixel_solid_angle(nside)
    t_pix = t_mission * n_det * omega_pix / OMEGA_SKY
    return t_pix


def time_per_arcmin_sq(t_mission=T_MISSION, n_det=1):
    """
    Time spent by a detector per arcmin square of the sky during a period t_mission and assuming an uniform scan.
    Default units : seconds
    """
    t_arcmin_sq = t_mission * n_det * OMEGA_ARCMIN_SQ / OMEGA_SKY
    return t_arcmin_sq


def time_per_solid_angle(t_mission=T_MISSION, n_det=1):
    """
    Time spent by a detector per unit solid angle of the sky during a period t_mission and assuming an uniform scan.
    Unit solid angle = 1 steradian
    Default units : seconds
    """
    t_unit_solid_angle = t_mission * n_det * 1.0 / OMEGA_SKY
    return t_unit_solid_angle


def sensitivity_to_noise_pix(det_sens, nside, t_mission=T_MISSION, n_det=N_DET_145):
    """
    Converts from the detector sensitivity to the noise in each pixel of a Healpix map with the given nside, assuming an uniform scan over the time of t_mission and by n_det detectors.
    det_sens is a single scalar value of unit uK.sqrt(s)
    Array of 3 values returned for T, Q, U maps respectively. The noise is sqrt(2) times more for Q and U
    Output unit : uK.arcmin
    """
    t_pix = time_per_pixel(nside, t_mission)
    noise_pix_T = det_sens / np.sqrt(t_pix) / np.sqrt(n_det)
    factor = np.array([1.0, np.sqrt(2), np.sqrt(2)])
    return noise_pix_T*factor


def sensitivity_to_noise_arcmin(det_sens, t_mission=T_MISSION, n_det=N_DET_145):
    """
    Converts from the detector sensitivity to the noise in each arcmin square a Healpix map, assuming an uniform scan over the time of t_mission and by n_det detectors.
    det_sens is a single scalar value of unit uK.sqrt(s)
    Array of 3 values returned for T, Q, U maps respectively. The noise is sqrt(2) times more for Q and U
    Output unit : uK.arcmin
    """
    t_arcmin_sq = time_per_arcmin_sq(t_mission)
    noise_arcmin_T = det_sens / np.sqrt(t_arcmin_sq) / np.sqrt(n_det)
    factor = np.array([1.0, np.sqrt(2), np.sqrt(2)])
    return noise_arcmin_T*factor


def sensitivity_to_noise_solid_angle(det_sens, t_mission=T_MISSION, n_det=N_DET_145):
    """
    Converts from the detector sensitivity to the noise in each unit solid angle of a Healpix map, assuming an uniform scan over the time of t_mission and by n_det detectors.
    det_sens is a single scalar value of unit uK.sqrt(s)
    Array of 3 values returned for T, Q, U maps respectively. The noise is sqrt(2) times more for Q and U
    Output unit : uK.arcmin
    """
    t_unit_solid_angle = time_per_solid_angle(t_mission)
    noise_solid_angle_T = det_sens / np.sqrt(t_unit_solid_angle) / np.sqrt(n_det)
    factor = np.array([1.0, np.sqrt(2), np.sqrt(2)])
    return noise_solid_angle_T*factor


def noise_arcmin_to_noise_pix(noise_arcmin, nside): 
    """
    Converts from noise in a arcmin square of the sky to noise in each Healpix pixel for the given nside.
    Input units : muK.arcmin
    Output units : muK
    """
    omega_pix = pixel_solid_angle(nside)
    noise_pix = noise_arcmin*np.sqrt(OMEGA_ARCMIN_SQ/omega_pix)
    return noise_pix


def noise_arcmin_to_noise_solid_angle(noise_arcmin):
    """
    Converts from noise in a arcmin square of the sky to noise in a unit solid angle of the sky.
    Input units : muK.arcmin
    Output units : muK
    """
    noise_solid_angle = noise_arcmin*np.sqrt(OMEGA_ARCMIN_SQ)
    return noise_solid_angle


def noise_arcmin_to_sensitivity(noise_arcmin, t_mission=T_MISSION, n_det=N_DET_145, comp="T"):
    """
    Converts from the noise in each arcmin square a Healpix map to detector sensitivity, assuming an uniform scan over the time of t_mission and by n_det detectors.
    noise_arcmin is a single scalar value of unit uK.arcmin. The noise in a pixel is sqrt(2) times more for Q and U
    det_sens is a scalar value of unit uK.sqrt(s)
    """
    t_arcmin_sq = time_per_arcmin_sq(t_mission)
    det_sens = noise_arcmin * np.sqrt(t_arcmin_sq) * np.sqrt(n_det)
    return det_sens


def sensitivity_to_noise_maps(det_sens, nside, t_mission=T_MISSION, n_det=N_DET_145, pol=True): 
    """
    Converts from the detector sensitivity to the noise in each pixel of a Healpix map with the given nside, assuming an uniform scan over the time of t_mission and by n_det detectors. It them makes a sky map with the calculated noise level.
    det_sens is a single scalar value of unit uK.sqrt(s)
    Array of 3 values returned for T, Q, U maps respectively if pol is True, otherwise only a T map is made. The noise is sqrt(2) times more for Q and U
    """
    npix = 12*nside**2
    if pol:
        noise_pix = sensitivity_to_noise_pix(det_sens, nside, t_mission, n_det)
        noise_map = np.empty((3, npix))
        for i in range(3):
            noise_map[i] = np.random.normal(loc=0.0, scale=noise_pix[i], size=npix)

    else:
        noise_pix = sensitivity_to_noise_pix(det_sens, nside, t_mission, n_det)[0]
        noise_map = np.random.normal(loc=0.0, scale=noise_pix, size=npix)

    return noise_map


def noise_arcmin_to_noise_maps(noise_arcmin, nside, pol=True):
    """
    Converts from the noise in arcmin square of the sky to the noise in each pixel of a Healpix map with the given nside. It them makes a sky map with the calculated noise level.
    noise_arcmin is a single or a set of 3 scalar values of unit uK.arcmin
    Array of 3 values returned for T, Q, U maps respectively if pol is True, otherwise only a T map is made. The noise is sqrt(2) times more for Q and U
    """
    npix = 12*nside**2
    if pol:
        sky = np.empty((3, npix))
        for i in range(3):
            noise_pix = noise_arcmin_to_noise_pix(noise_arcmin[i], nside)
            noise_map[i] = np.random.normal(loc=0.0, scale=noise_pix[i], size=npix)  
    else:
        noise_map = np.random.normal(loc=0.0, scale=noise_pix[0], size=npix)  
    return noise_map


def sensitivity_to_noise_mean_Cl(det_sens, beam_fwhm, lmax, t_mission=T_MISSION, n_det=N_DET_145, f_sky=1.0):
    """
    Given the sensitivity of the instrument, observation time and number of detectors, and assuming that the scan was uniform, tis gives us the theoretical mean of the noise spectra for TT, EE, and BB.
    """
    w_inv = sensitivity_to_noise_solid_angle(det_sens, t_mission, n_det)**2
    Bl_squared = hp.gauss_beam(fwhm=beam_fwhm, lmax=lmax, pol=True)**2
    noise_mean = np.empty((3, lmax+1))
    for i in range(3):
        noise_mean[i] = np.sqrt(1.0/f_sky)*w_inv[i]/Bl_squared[...,i]
        noise_mean[i][:2] = 0.0
    return noise_mean


def cosmic_variance_Cl(lmax, f_sky=1.0, spectra_fid=None):
    """
    This gives us the variance on the fiducial CMB power spectra due to cosmic variance alone.
    """
    if spectra_fid == None:
        spectra_fid = spectra_file
    Cl = np.load(spectra_fid)[..., :lmax+1]
    ell = np.arange(lmax + 1)
    variance = np.empty((3, lmax+1))
    for i in range(3):
        variance[i] = (2.0/(2.0*ell + 1)/f_sky)*Cl[i]**2
        variance[i][:2] = 0.0
    return variance


def sensitivity_to_noise_variance_Cl(det_sens, beam_fwhm, lmax, t_mission=T_MISSION, n_det=N_DET_145, f_sky=1.0):
    """
    Given the sensitivity of the instrument, observation time and number of detectors, and assuming that the scan was uniform, this gives us the theoretical variance of the noise spectra for TT, EE, and BB.
    """
    w_inv = sensitivity_to_noise_solid_angle(det_sens, t_mission, n_det)**2
    ell = np.arange(lmax + 1)
    Bl_squared = hp.gauss_beam(fwhm=beam_fwhm, lmax=lmax, pol=True)**2
    variance = np.empty((3, lmax+1))
    for i in range(3):
        variance[i] = (2.0/(2.0*ell + 1)/f_sky)*(w_inv[i]/Bl_squared[...,i])**2
        variance[i][:2] = 0.0
    return variance


def sensitivity_to_variance_Cl(det_sens, beam_fwhm, lmax, t_mission=T_MISSION, n_det=N_DET_145, f_sky=1.0, spectra_fid=None):
    """
    Given the sensitivity of the instrument, observation time and number of detectors, and assuming that the scan was uniform, this gives us the theoretical variance of the noise spectra and the cosmic variance for TT, EE, and BB.
    """
    if spectra_fid == None:
        spectra_fid = spectra_file
    Cl = np.load(spectra_fid)[..., :lmax+1]
    ell = np.arange(lmax + 1)
    w_inv = sensitivity_to_noise_solid_angle(det_sens, t_mission, n_det)**2
    Bl_squared = hp.gauss_beam(fwhm=beam_fwhm, lmax=lmax, pol=True)**2
    variance = np.empty((3, lmax+1))
    for i in range(3):
        variance[i] = (2.0/(2.0*ell + 1)/f_sky)*(Cl[i] + w_inv[i]/Bl_squared[...,i])**2
        variance[i][:2] = 0.0
    return variance
