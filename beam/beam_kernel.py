#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
from simulation.lib.plotting.my_imshow import new_imshow


def gaussian_2d(settings, mesh):
    factor = 2*np.sqrt(2*np.log(2))
    sigma_major = settings.fwhm_major/factor
    sigma_minor = settings.fwhm_minor/factor
    x0, y0 = settings.center
    theta = np.deg2rad(settings.tilt)
    x,y = mesh
    a = (np.cos(theta)**2)/(2*sigma_major**2) + (np.sin(theta)**2)/(2*sigma_minor**2)
    b = -1*np.sin(2*theta)/(4*sigma_major**2) + np.sin(2*theta)/(4*sigma_minor**2)
    c = (np.sin(theta)**2)/(2*sigma_major**2) + (np.cos(theta)**2)/(2*sigma_minor**2)
    beam_kernel = np.exp(-1*(a*(x - x0)**2 + 2*b*(x - x0)*(y - y0) + c*(y - y0)**2)) 
    if settings.normalise_beam:
        norm = 2*np.pi*sigma_major*sigma_minor
        beam_kernel /= norm
    return beam_kernel

def get_mesh(settings):
    factor = 2*np.sqrt(2*np.log(2))
    sigma_max = max(settings.fwhm_major, settings.fwhm_minor)/factor
    offset_max = max(np.abs(settings.center))
    size = settings.beam_cutoff*sigma_max + offset_max
    dd = settings.beam_resolution
    nxy = int(size/dd)
    x = np.linspace(-1.0*nxy*dd, nxy*dd, 2*nxy + 1)
    y = np.linspace(-1.0*nxy*dd, nxy*dd, 2*nxy + 1)
    return np.meshgrid(x,y), x 

def convolve_mesh(mesh, settings):
    pass

def display_beam_settings():
    if settings.do_pencil_beam:
        print "Pencil beam"
    else:
        print "Major axis(FWHM) : ", settings.fwhm_major, "arcmins"
        print "Minor axis(FWHM) : ", settings.fwhm_minor, "arcmins"
        print "Center : ", settings.center
        print "Tilt : ", settings.tilt, " degrees"
        print "Pixel size : ", settings.beam_resolution, "arcmins" 
        print "No of pixels per fwhm (minor-axis): ", 2*settings.fwhm_minor/settings.beam_resolution

def plot_beam():
    fig, ax = plt.subplots()
    im = new_imshow(ax, beam_kernel)
    fig.colorbar(im, ax=ax)
    plt.show()


if __name__=="__main__":
    from custom_settings import settings
    if settings.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        mesh, del_beta = get_mesh(settings)
        beam_kernel = gaussian_2d(settings, mesh)

    if settings.display_beam_settings:
        display_beam_settings()
    if settings.plot_beam:
        plot_beam()

def get_beam(settings=None):
    if settings is None:
        from custom_settings import settings
    if settings.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        mesh, del_beta = get_mesh(settings)
        beam_kernel = gaussian_2d(settings, mesh)

    return beam_kernel, del_beta
