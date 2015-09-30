#! /usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import pycamb

def make_spectra(settings = None):
    if settings is None:
        from settings_spectrum import settings
    if settings.cosmo_params['source'] is "pycamb_default":
        spectra = pycamb.camb(settings.lmax)
    elif settings.cosmo_params['source'] is "planck_latest":
        spectra = pycamb.camb(settings.lmax, settings.cosmo_params)
    else:
        spectra = pycamb.camb(settings.lmax, settings.cosmo_params)
        
    ell = np.arange(2, settings.lmax + 1)

    if settings.plot_spectra:
        plt.figure()
        make_subplot(1, 'TT', spectra[0], ell, '', r'$\mu K^2$')    
        make_subplot(2, 'EE', spectra[1], ell, '', '') 
        make_subplot(3, 'BB', spectra[2], ell, 'l', r'$\mu K^2$')    
        make_subplot(4, 'TE', spectra[3], ell, 'l', '')    
        plt.show()

    if settings.normalise_spectra is False:
        norm_fact = ell*(ell + 1)/(2*np.pi)
        spectra = spectra/norm_fact

    if settings.write_spectra:
        np.save(settings.spectra_file, spectra)
    
    if settings.return_spectra:
        return spectra

def make_subplot(plot_number, title, plot, ell, xlabel, ylabel):
    plt.subplot(2,2,plot_number)
    plt.plot(ell, plot, 'r-')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.tick_params(axis='x', labelsize=8)
    plt.tick_params(axis='y', labelsize=8)
    plt.annotate(title, xy=(0.85, 0.85), xycoords='axes fraction')

if __name__=="__main__":
    from settings_spectrum import settings
    make_spectra(settings)
