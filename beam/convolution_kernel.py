#! /usr/bin/env python 

import numpy as np
import healpy as hp
from scipy import ndimage
import matplotlib.pyplot as plt
import sys
from simulation.lib.plotting.my_imshow import new_imshow


def gaussian_2d(beam_params, bolo_params, mesh):
    factor = 2*np.sqrt(2*np.log(2))
    sigma_major = bolo_params.fwhm_major/factor
    sigma_minor = bolo_params.fwhm_minor/factor
    sigma = np.sqrt(sigma_major**2 - sigma_minor**2)
    theta = np.deg2rad(bolo_params.beam_angle)
    x,y = mesh
    beam_kernel = np.exp(-x**2/(2*sigma**2))
    dim = y[0].size
    beam_kernel[:dim/2] = 0
    beam_kernel[dim/2+1:] = 0
    beam_kernel[dim/2, dim/2] = 0
    #beam_kernel = ndimage.interpolation.rotate(beam_kernel, angle=30.0)
    normalisation_factor = np.sqrt(2*np.pi*sigma**2)
    beam_kernel /= normalisation_factor
    return beam_kernel


def check_normalisation(beam_params, bolo_params, beam_kernel):
    dx = beam_params.beam_resolution
    dy = beam_params.beam_resolution
    integral = np.sum(beam_kernel)*dx*dy
    print "The integral of the beam is :", integral
    print "Percentage difference with unity :", 100*(1-integral)


def get_mesh(beam_params, bolo_params):
    offset_max = max(abs(bolo_params.del_x), abs(bolo_params.del_y))              #arc-mins
    fwhm_major = bolo_params.fwhm_major
    fwhm_minor = bolo_params.fwhm_minor
    fwhm = np.sqrt(fwhm_major**2 - fwhm_minor**2)
    #fwhm_minor = bolo_params.fwhm
    #fwhm_major = fwhm_minor/np.sqrt(1 - bolo_params.ellipticity**2)
    size = beam_params.beam_cutoff*fwhm + offset_max             #arc-mins
    dd = beam_params.beam_resolution                                            #arc-mins
    nxy = int(size/dd/2)
    x = np.arange(-nxy, nxy+1)*dd
    y = -1*np.arange(-nxy, nxy+1)*dd
    return np.meshgrid(x,y), x 


def display_beam_settings(beam_params, bolo_params, mesh):
    if beam_params.do_pencil_beam:
        print "Pencil beam"
    else:
        factor = 2*np.sqrt(2*np.log(2))
        fwhm_major = bolo_params.fwhm_major
        fwhm_minor = bolo_params.fwhm_minor
        #fwhm_minor = bolo_params.fwhm
        #fwhm_major = fwhm_minor/np.sqrt(1 - bolo_params.ellipticity**2)
        print "Major axis(FWHM) :", fwhm_major, "arcmins" 
        print "Minor axis(FWHM) :", fwhm_minor, "arcmins"
        print "Center :", bolo_params.del_x, bolo_params.del_y
        print "Tilt :", bolo_params.beam_angle, "degrees"
        print "Pixel size :", beam_params.beam_resolution, "arcmins" 
        print "Kernel width in FWHM of beam:", beam_params.beam_cutoff
        print "# of pixels per FWHM (minor-axis) of beam :", fwhm_minor/beam_params.beam_resolution
        print "Expected # of pixels in kernel cross-section :", int(beam_params.beam_cutoff*fwhm_major/beam_params.beam_resolution/2)*2 + 1 
        print "Actual # of pixels in kernel cross-section :", mesh[0][0].size 

def plot_beam(beam_kernel):
    fig, ax = plt.subplots()
    im = new_imshow(ax, beam_kernel)
    fig.colorbar(im, ax=ax)
    plt.show()


if __name__=="__main__":

    from custom_params import beam_params
    from simulation.timestream_simulation.bolo_params.bolo_0001 import bolo_params

    if beam_params.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        mesh, del_beta = get_mesh(beam_params, bolo_params)
        beam_kernel = gaussian_2d(beam_params, bolo_params, mesh)

    if beam_params.check_normalisation:
        check_normalisation(beam_params, bolo_params, beam_kernel)
    if beam_params.display_beam_settings:
        display_beam_settings(beam_params, bolo_params, mesh)

    #beam_kernel/=np.max(beam_kernel)

    if beam_params.plot_beam:
        plot_beam(beam_kernel)


def get_beam(beam_params, bolo_params):
    if beam_params.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        mesh, del_beta = get_mesh(beam_params, bolo_params)
        beam_kernel = gaussian_2d(beam_params, bolo_params, mesh)
    return beam_kernel, del_beta
