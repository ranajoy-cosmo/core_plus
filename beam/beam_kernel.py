#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
from simulation.lib.plotting.my_imshow import new_imshow

class Beam:

    def __init__(self, settings=None):
        self.settings = settings

    def gaussian_2d(self):
        factor = 2*np.sqrt(2*np.log(2))
        sigma_major = self.settings.fwhm_major/factor
        sigma_minor = self.settings.fwhm_minor/factor
        x0, y0 = self.settings.center
        theta = np.deg2rad(settings.tilt)
        x,y = self.mesh
        a = (np.cos(theta)**2)/(2*sigma_major**2) + (np.sin(theta)**2)/(2*sigma_minor**2)
        b = -1*np.sin(2*theta)/(4*sigma_major**2) + np.sin(2*theta)/(4*sigma_minor**2)
        c = (np.sin(theta)**2)/(2*sigma_major**2) + (np.cos(theta)**2)/(2*sigma_minor**2)
        self.kernel = np.exp(-1*(a*(x - x0)**2 + 2*b*(x - x0)*(y - y0) + c*(y - y0)**2)) 
        if self.settings.normalise_beam:
            norm = 2*np.pi*sigma_major*sigma_minor
            self.kernel /= norm

    def get_mesh(self):
        factor = 2*np.sqrt(2*np.log(2))
        sigma_max = max(self.settings.fwhm_major, self.settings.fwhm_minor)/factor
        offset_max = max(np.abs(self.settings.center))
        size = self.settings.beam_cutoff*sigma_max + offset_max
        dd = self.settings.beam_resolution
        nxy = int(size/dd)
        x = np.linspace(-1.0*nxy*dd, nxy*dd, 2*nxy + 1)
        y = np.linspace(-1.0*nxy*dd, nxy*dd, 2*nxy + 1)
        self.del_beta = x 
        self.mesh = np.meshgrid(x,y)

    def display_beam_settings(self):
        if self.settings.do_pencil_beam:
            print "Pencil beam"
        else:
            print "Major axis(FWHM) : ", self.settings.fwhm_major, "arcmins"
            print "Minor axis(FWHM) : ", self.settings.fwhm_minor, "arcmins"
            print "Center : ", self.settings.center
            print "Tilt : ", self.settings.tilt, " degrees"
            print "Pixel size : ", self.settings.beam_resolution, "arcmins" 
            print "No, of pixels per fwhm (minor-axis): ", 2*self.settings.fwhm_minor/self.settings.beam_resolution

    def plot_beam(self):
        fig, ax = plt.subplots()
        im = new_imshow(ax, self.kernel)
        fig.colorbar(im, ax=ax)
        plt.show()


if __name__=="__main__":
    from local_settings import settings
    if settings.do_pencil_beam:
        beam = Beam(settings)
        beam.kernel = np.array([[1]])
    else:
        beam = Beam(settings)
        beam.get_mesh()
        beam.gaussian_2d()

    if settings.display_beam_settings:
        beam.display_beam_settings()
    if settings.plot_beam:
        beam.plot_beam()
else:
    from local_settings import settings
    if settings.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        beam = Beam(settings)
        beam.get_mesh()
        beam.gaussian_2d()
        beam_kernel = beam.kernel
        del_beta = beam.del_beta
