#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os
import sys
import simulation.power_spectrum.noise_power_spectra as nps

def plot_histogram(hist_number, num_bins, range, log=False, histtype='step', alpha=0.5, colours=['b','g','r','c','m','y']):
    plt.figure()
    weights = np.ones(12*nside**2) / (12*nside**2)
    plt.hist(0.25 * noise_rms**2 * pix_area * cov_map_1_1[hist_number], bins=num_bins, range=range, weights=weights, label="1 detector : 0.5 year", histtype=histtype, log=log, alpha=alpha, color=colours[0], zorder=10)
    plt.hist(0.25 * noise_rms**2 * pix_area * cov_map_1_2[hist_number], bins=num_bins, range=range, weights=weights, label="1 detector : 1 year", histtype=histtype, log=log, alpha=alpha, color=colours[1], zorder=8)
    plt.hist(0.25 * noise_rms**2 * pix_area * cov_map_1_3[hist_number], bins=num_bins, range=range, weights=weights, label="1 detector : 4 years", histtype=histtype, log=log, alpha=alpha, color=colours[2], zorder=6)
    plt.hist(0.25 * noise_rms**2 * pix_area * cov_map_2_1[hist_number], bins=num_bins, range=range, weights=weights, label="2 detectors : 0.5 year", histtype=histtype, log=log, alpha=alpha, color=colours[3], zorder=9)
    plt.hist(0.25 * noise_rms**2 * pix_area * cov_map_2_2[hist_number], bins=num_bins, range=range, weights=weights, label="2 detectors : 1 year", histtype=histtype, log=log, alpha=alpha, color=colours[4], zorder=7)
    plt.hist(0.25 * noise_rms**2 * pix_area * cov_map_2_3[hist_number], bins=num_bins, range=range, weights=weights, label="2 detectors : 4 years", histtype=histtype, log=log, alpha=alpha, color=colours[5], zorder=5)
    if hist_number == 0:
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 0.5*365.25*24*60*60.0, 1)[0]**2
        plt.axvline(x=avg_sens, color=colours[0])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 1*365.25*24*60*60.0, 1)[0]**2
        plt.axvline(x=avg_sens, color=colours[1])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 4*365.25*24*60*60.0, 1)[0]**2
        plt.axvline(x=avg_sens, color=colours[2])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 0.5*365.25*24*60*60.0, 2)[0]**2
        plt.axvline(x=avg_sens, color=colours[3])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 1*365.25*24*60*60.0, 2)[0]**2
        plt.axvline(x=avg_sens, color=colours[4])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 4*365.25*24*60*60.0, 2)[0]**2
        plt.axvline(x=avg_sens, color=colours[5])
    if hist_number == 3 or hist_number == 5:
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 0.5*365.25*24*60*60.0, 1)[1]**2
        plt.axvline(x=avg_sens, color=colours[0])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 1*365.25*24*60*60.0, 1)[1]**2
        plt.axvline(x=avg_sens, color=colours[1])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 4*365.25*24*60*60.0, 1)[1]**2
        plt.axvline(x=avg_sens, color=colours[2])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 0.5*365.25*24*60*60.0, 2)[1]**2
        plt.axvline(x=avg_sens, color=colours[3])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 1*365.25*24*60*60.0, 2)[1]**2
        plt.axvline(x=avg_sens, color=colours[4])
        avg_sens = nps.sensitivity_to_noise_arcmin(det_sens, 4*365.25*24*60*60.0, 2)[1]**2
        plt.axvline(x=avg_sens, color=colours[5])

    plt.xlabel("$[\mu K^2.arcmin^2]$", fontsize=15)
    plt.ylabel("$f_{sky}$", fontsize=15)
    plt.legend()
    plt.show()

if __name__=="__main__":
    map_dir = sys.argv[1]
    det_sens = float(sys.argv[2])
    sampling_rate = float(sys.argv[3])
    nside = float(sys.argv[4])

    out_dir = "/global/homes/b/banerji/simulation/output"

    cov_map_1_1 = hp.read_map(os.path.join(out_dir, map_dir, "1_det_half_year", "covariance_maps.fits"), field=(0,1,2,3,4,5))
    cov_map_1_2 = hp.read_map(os.path.join(out_dir, map_dir, "1_det_one_year", "covariance_maps.fits"), field=(0,1,2,3,4,5))
    cov_map_1_3 = hp.read_map(os.path.join(out_dir, map_dir, "1_det_four_years", "covariance_maps.fits"), field=(0,1,2,3,4,5))
    cov_map_2_1 = hp.read_map(os.path.join(out_dir, map_dir, "2_det_half_year", "covariance_maps.fits"), field=(0,1,2,3,4,5))
    cov_map_2_2 = hp.read_map(os.path.join(out_dir, map_dir, "2_det_one_year", "covariance_maps.fits"), field=(0,1,2,3,4,5))
    cov_map_2_3 = hp.read_map(os.path.join(out_dir, map_dir, "2_det_four_years", "covariance_maps.fits"), field=(0,1,2,3,4,5))

    noise_rms = det_sens * np.sqrt(sampling_rate)

    pix_area = hp.nside2resol(nside, arcmin=True)**2
