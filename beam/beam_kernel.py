#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from simulation.lib.plotting.my_imshow import new_imshow
from settings import settings

def gaussian_2d(x, y):
    factor = 2*np.sqrt(2*np.log(2))
    sigma_major = settings.fwhm_major/factor
    sigma_minor = settings.fwhm_minor/factor
    x0, y0 = settings.center
    theta = np.deg2rad(settings.tilt)
    a = (np.cos(theta)**2)/(2*sigma_major**2) + (np.sin(theta)**2)/(2*sigma_minor**2)
    b = -1*np.sin(2*theta)/(4*sigma_major**2) + np.sin(2*theta)/(4*sigma_minor**2)
    c = (np.sin(theta)**2)/(2*sigma_major**2) + (np.cos(theta)**2)/(2*sigma_minor**2)
    beam = np.exp(-1*(a*(x - x0)**2 + 2*b*(x - x0)*(y - y0) + c*(y - y0)**2)) 
    if settings.normalise_beam:
        norm = 2*np.pi*sigma_major*sigma_minor
        return (1/norm)*beam
    else:
        return beam

def get_mesh():
    factor = 2*np.sqrt(2*np.log(2))
    fwhm_max = np.max((settings.fwhm_major, settings.fwhm_minor))
    sigma = fwhm_max/factor
    dd = hp.nside2resol(settings.nside, arcmin = True)
    nxy = int(settings.beam_cutoff*sigma/dd)
    x = np.linspace(-1.0*nxy*dd, nxy*dd, 2*nxy + 1)
    y = np.linspace(-1.0*nxy*dd, nxy*dd, 2*nxy + 1)
    print x.size, y.size
    return np.meshgrid(x,y), y

def display_beam_settings():
    print "Major axis(FWHM) : ", settings.fwhm_major, "arcmins"
    print "Minor axis(FWHM) : ", settings.fwhm_minor, "arcmins"
    print "Center : ", settings.center
    print "Tilt : ", settings.tilt, " degrees"
    print "NSIDE : ", settings.nside
    print "Pixel size : ", hp.nside2resol(settings.nside, arcmin = True)
    print "No, of pixels per fwhm (minor-axis): ", 2*settings.fwhm_minor/hp.nside2resol(settings.nside, arcmin = True)

def plot_beam(beam):
    fig, ax = plt.subplots()
    im = new_imshow(ax, beam)
    fig.colorbar(im, ax=ax)
    plt.show()


if __name__=="__main__":
    (xx, yy), del_beta = get_mesh()
    beam_kernel = gaussian_2d(xx, yy)
    if settings.display_beam_settings:
        display_beam_settings()
    if settings.plot_beam:
        plot_beam(beam_kernel)

else:
    (xx, yy), del_beta = get_mesh()
    del_beta = np.deg2rad(del_beta/60.0)
    beam_kernel = gaussian_2d(xx, yy)
