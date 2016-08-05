#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_theoretical(lmax, plot_log=True, plot_list=["TT", "EE", "BB"]):
    spectra_folder = "/global/homes/b/banerji/simulation/spectra/"

    ell = np.arange(lmax+1)

    if plot_log:
        plotter = plt.loglog
    else:
        plotter = plt.plot

    spectra = np.load(spectra_folder + "r_0001/unlensed_cls.npy")[..., :lmax+1]

    if "TT" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[0]/2/np.pi, color='k')
        plt.annotate("TT", xy=(1.25, 950), size=12)

    if "EE" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[1]/2/np.pi, color='k')
        plt.annotate("EE", xy=(1.25, 0.039), size=12)

    if "BB" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[2]/2/np.pi, color='k')

    spectra = np.load(spectra_folder + "r_001/unlensed_cls.npy")[..., :lmax+1]

    if "BB" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[2]/2/np.pi, color='k')

    spectra = np.load(spectra_folder + "r_001/lensedtot_cls.npy")[..., :lmax+1]

    if "BB" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[2]/2/np.pi, color='k')

    spectra = np.load(spectra_folder + "r_0001/lensedtot_cls.npy")[..., :lmax+1]

    if "BB" in plot_list:
        plotter(ell, ell*(ell+1)*spectra[2]/2/np.pi, color='k')

        plt.annotate("BB", xy=(1.25, 7.5e-5), size=12)

        plt.annotate("r = 0.01", xy=(1.25, 5.0e-4), size=10)

        plt.annotate("r = 0.001", xy=(1.25, 1.1e-5), size=10)

    plt.show()


def make_decorations(ylim=[1e-11, 100], xlim=[1,5000]):
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.xlabel('l')
    plt.ylabel('$l(l+1)C_l/2\pi$ $[\mu K^2]$')
    plt.legend(loc="lower right", prop={'size':12})
