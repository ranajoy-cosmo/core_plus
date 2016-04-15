#! /usr/bin/env python 

import numpy as np
import healpy as hp
from scipy import ndimage
import matplotlib.pyplot as plt
import sys
import importlib
from simulation.lib.plotting.my_imshow import new_imshow


def gaussian_2d(beam_params, bolo_params, mesh):
    factor = 2*np.sqrt(2*np.log(2))
    #sigma_minor = bolo_params.fwhm/factor
    #sigma_major = (1+bolo_params.ellipticity)*sigma_minor 
    #sigma = np.sqrt(sigma_major**2 - sigma_minor**2)
    sigma = bolo_params.conv_fwhm/factor
    theta = np.deg2rad(bolo_params.beam_angle)
    x,y = mesh
    conv_kernel = np.exp(-y**2/(2*sigma**2))
    dim = y[0].size
    conv_kernel[...,:dim/2] = 0
    conv_kernel[...,dim/2+1:] = 0
    conv_kernel = ndimage.interpolation.rotate(conv_kernel, angle=bolo_params.beam_angle + bolo_params.pol_ang, reshape=False)
    conv_kernel[conv_kernel<0] = 0
    integral = np.sum(conv_kernel)*beam_params.beam_resolution
    #integral = np.sqrt(2*np.pi)*sigma
    return conv_kernel/integral


def check_normalisation(beam_params, bolo_params, conv_kernel):
    integral = np.sum(conv_kernel)*beam_params.beam_resolution
    print "The integral of the beam is :", integral
    print "Percentage difference with unity :", 100*(1-integral)


def get_mesh(beam_params, bolo_params):
    #fwhm_minor = bolo_params.fwhm
    #fwhm_major = (1+bolo_params.ellipticity)*fwhm_minor 
    #fwhm = np.sqrt(fwhm_major**2 - fwhm_minor**2)
    fwhm = bolo_params.conv_fwhm
    size = beam_params.beam_cutoff*fwhm                                      #arc-mins
    dd = beam_params.beam_resolution                                            #arc-mins
    n = int(size/dd/2)
    x = np.arange(-n, n+1)*dd
    y = -1*np.arange(-n, n+1)*dd
    return np.meshgrid(x,y), x 


def display_beam_settings(beam_params, bolo_params, beam_kernel):
    fwhm_minor = bolo_params.fwhm
    fwhm_major = np.sqrt(fwhm_minor**2 + bolo_params.conv_fwhm**2)
    fwhm_conv = bolo_params.conv_fwhm
    ellipticity = 100*2*(fwhm_major - fwhm_minor)/(fwhm_major + fwhm_minor)
    print "Major axis(FWHM) :", fwhm_major, "arcmins" 
    print "Minor axis(FWHM) :", fwhm_minor, "arcmins"
    print "Convolution function FWHM :", fwhm_conv, "arcmins"
    print "Ellipticity :", ellipticity
    print "Tilt :", bolo_params.beam_angle, "degrees"
    print "Pixel size :", beam_params.beam_resolution, "arcmins" 
    print "Actual # of pixels in kernel cross-section :", beam_kernel[0].size 

def plot_beam(beam_kernel, beam_params):
    fig, ax = plt.subplots()
    n = beam_kernel[0].size/2
    extent = np.arange(-n, n+1)*beam_params.beam_resolution
    im = new_imshow(ax, beam_kernel, x=extent, y=extent, interpolation="nearest")
    fig.colorbar(im, ax=ax)
    plt.show()


if __name__=="__main__":

    from custom_params import beam_params
    bolo_name = sys.argv[1]
    bolo_params = importlib.import_module("simulation.timestream_simulation.bolo_params." + bolo_name).bolo_params

    if bolo_params.conv_fwhm == 0.0:
        conv_kernel = np.array([[1]])/beam_params.beam_resolution
        del_beta = np.array([0])
        beam_params.plot_beam=False
    else:
        mesh, del_beta = get_mesh(beam_params, bolo_params)
        conv_kernel = gaussian_2d(beam_params, bolo_params, mesh)

    if beam_params.display_beam_settings:
        display_beam_settings(beam_params, bolo_params, conv_kernel)
    if beam_params.check_normalisation:
        check_normalisation(beam_params, bolo_params, conv_kernel)
    if beam_params.plot_beam:
        plot_beam(conv_kernel, beam_params)


def get_beam(beam_params, bolo_params):
    mesh, del_beta = get_mesh(beam_params, bolo_params)
    conv_kernel = gaussian_2d(beam_params, bolo_params, mesh)
    return conv_kernel, del_beta
