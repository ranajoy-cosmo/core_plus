#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Path names
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

tag = "b_001_hybrid"#"b_001"
output_folder = "/global/u1/b/banerji/leap/apps/transfer_function/xpure/spectrum/"
input_spectra_folder = "/global/u1/b/banerji/leap/apps/simulated_timestreams/bolo/maps_and_spectra/spectra/"
output_plot_folder = "/global/u1/b/banerji/leap/apps/transfer_function/xpure/plots/"
binary_mask_folder = "/global/u1/b/banerji/leap/apps/simulated_timestreams/bolo/maps_and_spectra/masks/"
spectra_name = "spectra001.npy"#"spectra001.npy"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

LMAX = 1500
NSIDE = 1024
FWHM_arcmin = 8.0
FWHM = np.deg2rad(FWHM_arcmin/60.0)
SIGMA = FWHM/2.35482
ell = np.arange(LMAX + 1)
norm = 2*np.pi
beam_attenuation = np.exp(-ell*(ell + 1)*SIGMA**2)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Loading input spectra and masks
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

cl = hp.mrdfits(output_folder + "cls/" + tag + "/cellpure_sky_map_mask1_0_0.fits")
ps = hp.mrdfits(output_folder + "cls/" + tag + "/pseudopure_sky_map_mask1_0_0.fits")
input_spectra = np.load(input_spectra_folder + spectra_name)[...,:LMAX + 1]
binary_mask = hp.read_map(binary_mask_folder + "mask_ebex.fits")
apodised_mask = hp.read_map(output_folder + "mask/" + tag + "/apodized_mask_I_" + tag + ".fits")


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Calculate the normalisation due to partial sky and other factors
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

sky_frac = float(np.sum(binary_mask))/hp.nside2npix(NSIDE)
norm_apodised = float(np.sum(apodised_mask**2))/hp.nside2npix(NSIDE)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Plotting sub routine
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def make_subplot(plot_number, title, p1, p2, p3, xlabel, ylabel):
    plt.subplot(2,2,plot_number)
    plt.plot(ell[50:], p1[50:LMAX + 1]*ell[50:]*(ell[50:] + 1)/norm, 'g.', markersize=2)
    plt.plot(ell[50:], p2[50:LMAX + 1]*ell[50:]*(ell[50:] + 1)/norm, 'r-')
    if p3 is not None:
        plt.plot(ell[50:], p3[50:LMAX + 1]*ell[50:]*(ell[50:] + 1)/norm, 'y-')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel) 
    plt.tick_params(axis='x', labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    plt.annotate(title, xy=(0.85, 0.85), xycoords='axes fraction')

def make_subplot_semilog(plot_number, title, p1, p2, p3, xlabel, ylabel):
    plt.subplot(2,2,plot_number)
    plt.semilogy(ell, p1[:LMAX + 1]*ell*(ell + 1)/norm, 'g.', markersize=2)
    plt.semilogy(ell, p2[:LMAX + 1]*ell*(ell + 1)/norm, 'r-')
    if p3 is not None:
        plt.semilogy(ell, p3[:LMAX + 1]*ell*(ell + 1)/norm, 'y-')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel) 
    plt.tick_params(axis='x', labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    plt.annotate(title, xy=(0.85, 0.85), xycoords='axes fraction')

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#Making the xpure plots
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

expected_beamed_spectra = input_spectra*beam_attenuation

plt.figure(1)

make_subplot(1, 'TT', ps[0]/norm_apodised, input_spectra[0], expected_beamed_spectra[0], " ", r'$\mu K^2$')
plt.plot(cl[0][2:], cl[1][2:])
make_subplot(2, 'EE', ps[1]/norm_apodised, input_spectra[1], expected_beamed_spectra[1], " ", " ")
plt.plot(cl[0][2:], cl[2][2:])
make_subplot(3, 'BB', ps[2]/norm_apodised, input_spectra[2], expected_beamed_spectra[2], "l", r'$\mu K^2$')
plt.plot(cl[0][2:], cl[3][2:])
make_subplot(4, 'TE', ps[3]/norm_apodised, input_spectra[3], expected_beamed_spectra[3], "l", " ")
plt.plot(cl[0][2:], cl[4][2:])

plt.savefig(output_plot_folder + "xpure_output_" + tag + ".png")


plt.figure(2)

make_subplot_semilog(1, 'TT', ps[0]/norm_apodised, input_spectra[0], expected_beamed_spectra[0], " ", r'$\mu K^2$')
plt.semilogy(cl[0][:], cl[1][:])
make_subplot_semilog(2, 'EE', ps[1]/norm_apodised, input_spectra[1], expected_beamed_spectra[1], " ", " ")
plt.semilogy(cl[0][:], cl[2][:])
make_subplot_semilog(3, 'BB', ps[2]/norm_apodised, input_spectra[2], expected_beamed_spectra[2], "l", r'$\mu K^2$')
plt.semilogy(cl[0][:], cl[3][:])
make_subplot_semilog(4, 'TE', ps[3]/norm_apodised, input_spectra[3], expected_beamed_spectra[3], "l", " ")
plt.semilogy(cl[0][:], cl[4][:])

plt.savefig(output_plot_folder + "xpure_output_" + tag + "_semilog.png")


plt.figure(3)

plt.plot(cl[0][:25], cl[3][:25], 'b-')
plt.plot(ell[10:500], ell[10:500]*(ell[10:500] + 1)*ps[2][10:500]/norm/norm_apodised, 'g.', markersize=2)
plt.plot(ell[:500], ell[:500]*(ell[:500] + 1)*expected_beamed_spectra[2][:500]/norm, 'r-')

plt.savefig(output_plot_folder + "bb_comparison_" + tag + ".png")


plt.figure(4)

plt.semilogy(cl[0][:25], cl[3][:25], 'b-')
plt.semilogy(ell[10:500], ell[10:500]*(ell[10:500] + 1)*ps[2][10:500]/norm/norm_apodised, 'g.', markersize=2)
plt.semilogy(ell[:500], ell[:500]*(ell[:500] + 1)*expected_beamed_spectra[2][:500]/norm, 'r-')

plt.savefig(output_plot_folder + "bb_comparison_" + tag + "_semilog.png")

