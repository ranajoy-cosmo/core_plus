#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import sys
from simulation.lib.plotting.my_imshow import new_imshow
from simulation.lib.geometry.conversions import *


def get_beam_pixel_weights(dim, settings):
    del_theta = settings.beam_resolution/60.0
    del_phi = np.full(dim, settings.beam_resolution)/60.0
    thetas = np.arange(-dim/2+1, dim/2+1)*del_theta + settings.scan_radius
    del_phi = del_phi*np.sin(np.radians(thetas))/np.sin(np.radians(settings.scan_radius))
    
    central_pixel_area = del_theta*del_phi[dim/2]
    beam_weight = (del_theta*del_phi*np.ones(dim**2).reshape((dim,dim))/central_pixel_area).T

    return beam_weight

def gaussian_angular(settings, mesh):
    factor = 2*np.sqrt(2*np.log(2))
    sigma = settings.fwhm_major/factor
    beam_kernel = np.exp(-mesh**2/(2*sigma**2))
    dim = mesh[0].size
    beam_weights = get_beam_pixel_weights(dim, settings)
    return beam_kernel*beam_weights, beam_weights

def get_mesh(settings):
    size = settings.beam_cutoff*settings.fwhm_major                   #Half-width of the beam kernel in arcmin
    dd = settings.beam_resolution                       #Beam pixel size in arcmin
    n = int(size/dd/2)                                  #No. of pixels per half width

    lat = settings.scan_radius + np.arange(-n, n+1)*am2deg(dd)
    lon = np.arange(-n, n+1)*am2deg(dd)/np.sin(deg2rad(settings.scan_radius))

    llon, llat = np.meshgrid(lon, lat)

    return  (llon, llat), lat

def get_ang_distance(mesh):
    lon, lat = mesh
    dim = lon[0].size
    ci = dim/2                 #Central index

    centre = lon[ci, ci], lat[ci, ci]
    ang_dist = np.empty((dim,dim))

    for i in range(dim):
        for j in range(dim):
            point = lon[i, j], lat[i, j]
            ang_dist[i][j] = np.degrees(hp.rotator.angdist(point, centre, lonlat=True))*60.0

    return ang_dist

def display_beam_settings(settings, mesh):
    if settings.do_pencil_beam:
        print "Pencil beam"
    else:
        print "Major axis(FWHM) :", settings.fwhm_major, "arcmins"
        print "Minor axis(FWHM) :", settings.fwhm_minor, "arcmins"
        print "Center :", settings.center
        print "Tilt :", settings.tilt, "degrees"
        print "Pixel size :", settings.beam_resolution, "arcmins" 
        print "Kernel width in FWHM of beam:", settings.beam_cutoff 
        print "# of pixels per FWHM (minor-axis) of beam :", settings.fwhm_minor/settings.beam_resolution
        print "Expected # of pixels in kernel cross-section :", int(settings.beam_cutoff*settings.fwhm_major/settings.beam_resolution/2)*2 + 1 
        print "Actual # of pixels in kernel cross-section :", mesh[0][0].size 

def get_hitmap(mesh, beam_kernel):
    hitmap = np.zeros(12*settings.nside**2)
    pix = hp.ang2pix(settings.nside, np.radians(mesh[1].flatten()), np.radians(mesh[0].flatten()))
    hitmap[pix] = 1#np.arange(pix.size + 1)

    beam_healpix = np.zeros(12*settings.nside**2)
    beam_healpix[pix] = beam_kernel.flatten()

    return hitmap, beam_healpix


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
        ang_dist_mesh = get_ang_distance(mesh)
        beam_kernel, beam_weights = gaussian_angular(settings, ang_dist_mesh)
    beam_kernel/=np.max(beam_kernel)

    hitmap, beam_healpix = get_hitmap(mesh, beam_kernel)

    if settings.display_beam_settings:
        display_beam_settings(settings, mesh)

    if settings.plot_beam:
        plot_beam()

def get_beam(beam_params, bolo_params):
    if beam_params.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        mesh, del_beta = get_mesh(beam_params, bolo_params)
        ang_dist_mesh = get_ang_distance(mesh)
        beam_kernel, beam_weights = gaussian_angular(settings, ang_dist_mesh)
    beam_kernel/=np.max(beam_kernel)

    return beam_kernel, del_beta
