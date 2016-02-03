#! /usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import sys
import simulation.lib.quaternion.quaternion as qt
from simulation.lib.plotting.my_imshow import new_imshow
<<<<<<< HEAD:beam/beam_kernel_analytic_s2.py
from custom_settings import settings
from simulation.lib.geometry.conversions import *
=======
from simulation.lib.geometry.conversions import *

>>>>>>> 1dd1082e90be7068d7d9518d1c0d3ac20c426bbc:beam/beam_kernel_s2_v2.py

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
<<<<<<< HEAD:beam/beam_kernel_analytic_s2.py
    size = settings.beam_cutoff*settings.fwhm           #(arcmin)Full-width of the beam kernel in arcmin
    dd = settings.beam_resolution                       #Beam pixel size in arcmin
    n = int(size/dd)/2                                  #No. of pixels per half width

    lat = settings.scan_radius + np.arange(-n, n+1)*am2deg(dd)
    lon = np.arange(-n, n+1)*am2deg(dd)/np.sin(deg2rad(settings.scan_radius))
=======
    size = settings.beam_cutoff*settings.fwhm_major                   #Half-width of the beam kernel in arcmin
    dd = settings.beam_resolution                       #Beam pixel size in arcmin
    n = int(size/dd/2)                                  #No. of pixels per half width

    rot_axis = np.array([np.cos(deg2rad(settings.scan_radius), 0.0, np.sin(deg2rad(settings.scan_radius)))])
>>>>>>> 1dd1082e90be7068d7d9518d1c0d3ac20c426bbc:beam/beam_kernel_s2_v2.py

    angles = np.arange(-n, n+1)*am2deg(dd)
    mesh = np.full((angles.size, 2), 0)
    
    q = qt.make_quaternion(angles, axis, degree=True)

    for lat in angles:
        lon = qt.transform(q, np.array([])

<<<<<<< HEAD:beam/beam_kernel_analytic_s2.py
=======
    lon = np.arange(-n, n+1)*am2deg(dd)/np.sin(deg2rad(settings.scan_radius))

    llon, llat = np.meshgrid(lon, lat)

>>>>>>> 1dd1082e90be7068d7d9518d1c0d3ac20c426bbc:beam/beam_kernel_s2_v2.py
    return  (llon, llat), lat

def get_ang_distance(mesh):
    lon, lat = mesh
    dim = lon[0].size
    ci = dim/2                 #Central index

    centre = np.array([lon[ci, ci], lat[ci, ci]])
    #ang_dist = np.empty((dim*dim))

    #for i in range(dim):
    #    for j in range(dim):
    #        point = lon[i, j], lat[i, j]
    #        ang_dist[i][j] = np.degrees(hp.rotator.angdist(point, centre, lonlat=True))*60.0

    #ang_dist = hp.rotator.angdist(np.dstack((mesh[0], mesh[1])), centre, lonlat=True)
    ang_dist = hp.rotator.angdist(zip(deg2rad*lon.flatten(), deg2rad*lat.flatten()), deg2rad*centre, lonlat=True)
    ang_dist = rad2am(ang_dist).reshape((dim,dim))

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
<<<<<<< HEAD:beam/beam_kernel_analytic_s2.py
        print "Total # of pixels in kernel cross-section :", int(2*settings.beam_cutoff*settings.fwhm_major/factor/settings.beam_resolution) 

def get_hitmap(mesh, beam_kernel):
    hitmap = np.zeros(12*settings.nside**2)
    pix = hp.ang2pix(settings.nside, np.radians(mesh[1].flatten()), np.radians(mesh[0].flatten()))
    hitmap[pix] = 1#np.arange(pix.size + 1)

    beam_healpix = np.zeros(12*settings.nside**2)
    beam_healpix[pix] = beam_kernel.flatten()

    return hitmap, beam_healpix

def get_beam_pixel_weights(dim):
    del_theta = settings.beam_resolution/60.0
    del_phi = np.full(dim, settings.beam_resolution)/60.0
    thetas = np.arange(-dim/2+1, dim/2+1)*del_theta + settings.scan_radius
    del_phi = del_phi*np.sin(np.radians(thetas))/np.sin(np.radians(settings.scan_radius))
    
    central_pixel_area = del_theta*del_phi[dim/2]
    beam_weight = (del_theta*del_phi*np.ones(dim**2).reshape((dim,dim))/central_pixel_area).T

    return beam_weight

def check_integration(map1, map2, mesh, beam_kernel):
    dim = mesh[0][0].size
    centre = mesh[0][dim/2][dim/2], mesh[1][dim/2][dim/2]
    integral = 0

    beam_weight = get_beam_pixel_weights(dim)
    convolution = hp.get_interp_val(map1, np.radians(mesh[1].flatten()), np.radians(mesh[0].flatten()))*beam_kernel.flatten()*beam_weight.flatten()
    integral = np.sum(convolution)/np.sum(beam_kernel*beam_weight)

    output_nside = hp.get_nside(map2)
    input_nside = hp.get_nside(map1)
    input_pixel = hp.ang2pix(input_nside, np.radians(centre[1]), np.radians(centre[0]))
    output_pixel = hp.ang2pix(output_nside, np.radians(centre[1]), np.radians(centre[0]))
    print "Value before convolution :", map1[input_pixel]
    print "Integral using beam :", integral
    print "Value at pixel :", map2[output_pixel]
=======
        print "Expected # of pixels in kernel cross-section :", int(settings.beam_cutoff*settings.fwhm_major/settings.beam_resolution/2)*2 + 1 
        print "Actual # of pixels in kernel cross-section :", mesh[0][0].size 
>>>>>>> 1dd1082e90be7068d7d9518d1c0d3ac20c426bbc:beam/beam_kernel_s2_v2.py


def plot_beam():
    fig, ax = plt.subplots()
    im = new_imshow(ax, beam_kernel)
    fig.colorbar(im, ax=ax)
    plt.show()


if __name__=="__main__":
<<<<<<< HEAD:beam/beam_kernel_analytic_s2.py

    mesh, del_beta = get_mesh(settings)
    ang_dist_mesh = get_ang_distance(mesh)
    """
    beam_kernel = gaussian_angular(settings, ang_dist_mesh)
=======
    from custom_settings import settings
    if settings.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        mesh, del_beta = get_mesh(settings)
        ang_dist_mesh = get_ang_distance(mesh)
        beam_kernel, beam_weights = gaussian_angular(settings, ang_dist_mesh)
>>>>>>> 1dd1082e90be7068d7d9518d1c0d3ac20c426bbc:beam/beam_kernel_s2_v2.py
    beam_kernel/=np.max(beam_kernel)

    if settings.display_beam_settings:
        display_beam_settings(settings, mesh)

    if settings.plot_beam:
        plot_beam()
<<<<<<< HEAD:beam/beam_kernel_analytic_s2.py
        pass
    """
=======

def get_beam(settings=None):
    if settings is None:
        from custom_settings import settings
    if settings.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        mesh, del_beta = get_mesh(settings)
        ang_dist_mesh = get_ang_distance(mesh)
        beam_kernel, beam_weights = gaussian_angular(settings, ang_dist_mesh)
    beam_kernel/=np.max(beam_kernel)

    return beam_kernel, del_beta
>>>>>>> 1dd1082e90be7068d7d9518d1c0d3ac20c426bbc:beam/beam_kernel_s2_v2.py
