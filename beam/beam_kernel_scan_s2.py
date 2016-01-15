#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import sys, importlib
from simulation.lib.plotting.my_imshow import new_imshow
bolo_params = importlib.import_module("simulation.bolo.bolo_params.0001").bolo_params
from simulation.bolo.custom_settings import pointing_params 
from custom_settings import settings

def gaussian_angular(settings, mesh):
    factor = 2*np.sqrt(2*np.log(2))
    sigma = settings.fwhm_major/factor
    beam_kernel = np.exp(-mesh**2/(2*sigma**2))
    return beam_kernel

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

def generate_pointing(del_beta, nxy):
    beta = np.radians(bolo_params.beta + del_beta/60.0)
    u_init = np.array([np.cos(beta), 0.0, np.sin(beta)])
    w_prec = 2*np.pi/pointing_params.t_prec
    w_spin = 2*np.pi/pointing_params.t_spin
    t_steps = 0.001*(pointing_params.t_sampling/pointing_params.oversampling_rate)*np.arange(-nxy, nxy+1)
    R = po.Rotation3dOperator("XY'X''", -1.0*w_prec*t_steps, -1.0*np.full(2*nxy+1, np.radians(bolo_params.alpha)), -w_spin*t_steps)
    v = R*u_init
    #print v
    #print v.shape
    return v


def get_mesh(settings):
    factor = 2*np.sqrt(2*np.log(2))                     #Factor for conversion between FWHM and sigma
    sigma = settings.fwhm_major/factor                  #Corresponding sigma in arcmin
    size = settings.beam_cutoff*sigma                   #Half-width of the beam kernel in arcmin
    dd = settings.beam_resolution                       #Beam pixel size in arcmin
    nxy = int(size/dd)                                  #No. of pixels per half width

    if pointing_params.mode==1: 
        pointing_params.t_spin = 360.0*60.0*np.sin(bolo_params.beta)*pointing_params.t_sampling/1000.0/pointing_params.theta_co
        pointing_params.t_prec = 360.0*60.0*np.sin(bolo_params.alpha)*pointing_params.t_spin/pointing_params.theta_cross
    if pointing_params.mode==2:
        pointing_params.theta_cross = 360.0*60.0*np.sin(bolo_params.alpha)*pointing_params.t_spin/pointing_params.t_prec
        pointing_params.theta_co = 360*60*np.sin(bolo_params.beta)*pointing_params.t_sampling/1000.0/pointing_params.t_spin

    del_beta = np.arange(-nxy, nxy+1)*dd

    hitmap = np.zeros(12*settings.nside**2)

    llon = np.empty(0)
    llat = np.empty(0)
    for d_beta in del_beta:
        v = generate_pointing(d_beta, nxy)
        hitmap[hp.vec2pix(settings.nside, v[...,0], v[...,1], v[...,2])] = 1
        lat, lon = hp.vec2ang(v)
        llon = np.append(llon, lon)
        llat = np.append(llat, lat)
    llon = llon.reshape((23,23))
    llat = llat.reshape((23,23))
    return (llon, llat), del_beta, hitmap

def display_beam_settings():
    if settings.do_pencil_beam:
        print "Pencil beam"
    else:
        factor = 2*np.sqrt(2*np.log(2))
        print "Major axis(FWHM) :", settings.fwhm_major, "arcmins"
        print "Minor axis(FWHM) :", settings.fwhm_minor, "arcmins"
        print "Center :", settings.center
        print "Tilt :", settings.tilt, "degrees"
        print "Pixel size :", settings.beam_resolution, "arcmins" 
        print "Kernel width in FWHM of beam:", 2*settings.beam_cutoff/factor 
        print "# of pixels per FWHM (minor-axis) of beam :", settings.fwhm_minor/settings.beam_resolution
        print "Total # of pixels in kernel cross-section :", int(2*settings.beam_cutoff*settings.fwhm_major/factor/settings.beam_resolution) 

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
        mesh, del_beta, hitmap = get_mesh(settings)
        ang_dist_mesh = get_ang_distance(mesh)
        beam_kernel = gaussian_angular(settings, ang_dist_mesh)

    if settings.check_normalisation:
        #check_normalisation(settings, beam_kernel)
        pass
    if settings.display_beam_settings:
        display_beam_settings()
        pass
    beam_kernel/=np.max(beam_kernel)

    if settings.plot_beam:
        plot_beam()
        pass

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
