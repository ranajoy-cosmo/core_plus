#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from simulation.params.custom_params import global_paths

spectra_folder = os.path.join(global_paths.base_dir, "spectra")

def plot_theoretical(lmax, plot_log=True, plot_list=["TT", "EE", "BB"], lensed=True):

    ell = np.arange(lmax+1)

    if plot_log:
        plotter = plt.loglog
    else:
        plotter = plt.plot

    #Plotting the TT and EE spectra
    
    if lensed:
        spectra = np.load(os.path.join(spectra_folder, "r_0001/lensedtot_cls.npy"))[..., :lmax+1]
    else:
        spectra = np.load(os.path.join(spectra_folder, "r_0001/unlensed_cls.npy"))[..., :lmax+1]

    if "TT" in plot_list:
        plotter(ell, np.sqrt(ell*(ell+1)*spectra[0]/2/np.pi), color='k')
        plt.annotate("TT", xy=(1.25, 31), size=12)

    if "EE" in plot_list:
        plotter(ell, np.sqrt(ell*(ell+1)*spectra[1]/2/np.pi), color='k')
        plt.annotate("EE", xy=(1.25, 0.2), size=12)


    if "BB" in plot_list:
        spectra = np.load(os.path.join(spectra_folder, "r_0001/unlensed_cls.npy"))[..., :lmax+1]
        plotter(ell, np.sqrt(ell*(ell+1)*spectra[2]/2/np.pi), color='k')

        spectra = np.load(os.path.join(spectra_folder, "r_0001/lensedtot_cls.npy"))[..., :lmax+1]
        plotter(ell, np.sqrt(ell*(ell+1)*spectra[2]/2/np.pi), color='k')

        spectra = np.load(os.path.join(spectra_folder, "r_001/unlensed_cls.npy"))[..., :lmax+1]
        plotter(ell, np.sqrt(ell*(ell+1)*spectra[2]/2/np.pi), color='k')

        spectra = np.load(os.path.join(spectra_folder, "r_001/lensedtot_cls.npy"))[..., :lmax+1]
        plotter(ell, np.sqrt(ell*(ell+1)*spectra[2]/2/np.pi), color='k')

        plt.annotate("BB", xy=(1.25, 0.008), size=12)

        plt.annotate("r = 0.01", xy=(1.25, 0.019), size=10)

        plt.annotate("r = 0.001", xy=(1.25, 0.0035), size=10)

    plt.show()


def make_decorations(ylim=[1e-4, 100], xlim=[1,3000]):
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.xlabel('l')
    plt.ylabel('$\sqrt{l(l+1)C_l/2\pi}$ $[\mu K]$')
    plt.legend(loc="lower right", prop={'size':12})


def plot_spectra(cl, lmax=None, label=None, plot_log=True):
    if plot_log:
        plotter = plt.loglog
    else:
        plotter = plt.plot

    if lmax==None:
        lmax = cl.size - 1
    ell = np.arange(lmax+1)
    if label:
        plotter(ell, np.sqrt(ell*(ell+1)*cl[:lmax+1]/2/np.pi), label=label)
    else:
        plotter(ell, np.sqrt(ell*(ell+1)*cl[:lmax+1]/2/np.pi))
