#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_theoretical(lmax):
    spectra_folder = "/global/homes/b/banerji/simulation/spectra/"
    spectra = np.load(spectra_folder + "r_0001/unlensed_cls.npy")[..., :lmax+1]

    ell = np.arange(lmax+1)

    plt.loglog(ell, ell*(ell+1)*spectra[1]/2/np.pi, color='k')

    plt.loglog(ell, ell*(ell+1)*spectra[2]/2/np.pi, color='k')

    spectra = np.load(spectra_folder + "r_001/unlensed_cls.npy")[..., :lmax+1]

    plt.loglog(ell, ell*(ell+1)*spectra[2]/2/np.pi, color='k')

    spectra = np.load(spectra_folder + "r_001/lensedtot_cls.npy")[..., :lmax+1]

    plt.loglog(ell, ell*(ell+1)*spectra[2]/2/np.pi, color='k')

    spectra = np.load(spectra_folder + "r_0001/lensedtot_cls.npy")[..., :lmax+1]

    plt.loglog(ell, ell*(ell+1)*spectra[2]/2/np.pi, color='k')

    plt.annotate("EE", xy=(1.25, 0.039), size=12)

    plt.annotate("BB", xy=(1.25, 7.5e-5), size=12)

    plt.annotate("r = 0.01", xy=(1.25, 5.0e-4), size=10)

    plt.annotate("r = 0.001", xy=(1.25, 1.1e-5), size=10)

    plt.show()


def make_decorations():
    plt.ylim([1e-11,100])
    plt.xlim([1,5000])
    plt.xlabel('l')
    plt.ylabel('$l(l+1)C_l/2\pi$ $[\mu K^2]$')
    plt.legend(loc="lower right", prop={'size':12})
