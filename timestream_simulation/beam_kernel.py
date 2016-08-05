#! /usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.misc import imresize
import sys
import os
import importlib
from simulation.lib.plotting.my_imshow import new_imshow

class Beam():
    def __init__(self, config, bolo_config):
        self.config = config
        self.config.__dict__.update(bolo_config.__dict__)
        self.get_beam()


    def get_beam(self):
        if self.config.simulate_beam:
            if self.config.do_pencil_beam:
                self.beam_kernel = np.array([[1.0]])
                self.del_beta = np.array([0.0])
            else:
                mesh = self.get_mesh()
                self.gaussian_2d(mesh)
        else:
            self.read_mark_beam_map()


    def read_mark_beam_map(self):
        beam_kernel = np.load(self.config.beam_file)
        #orig_dim = 181                              #pixels
        #orig_res = 0.29796916162354475              #arc-mins 
        self.config.fwhm_major = 5.79
        self.config.fwhm_minor = 5.79
        mark_resolution = 0.2979691616235447
        num_pix = self.config.beam_cutoff * self.config.fwhm_major / mark_resolution
        n = int(num_pix/2)
        beam_kernel = beam_kernel[..., 91-n:91+n+1, 91-n:91+n+1]
        factor = self.config.beam_resolution / mark_resolution
        print factor
        print 2*n + 1
        new_num_pix = int((2*n+1)/factor)
        print new_num_pix
        self.beam_kernel = np.empty((4, new_num_pix, new_num_pix))
        for i in range(4):
            self.beam_kernel[i] = imresize(beam_kernel[i], (new_num_pix, new_num_pix))
            self.beam_kernel[i] /= (np.sum(self.beam_kernel[i])*self.config.beam_resolution**2)


    def gaussian_2d(self, mesh):
        factor = 2*np.sqrt(2*np.log(2))                                     #FWHM <--> sigma

        sigma_major = self.config.fwhm_major/factor
        sigma_minor = self.config.fwhm_minor/factor
        theta = np.deg2rad(self.config.beam_angle + self.config.pol_phase_ini)
        x0, y0 = self.config.offset_x/60.0, self.config.offset_y/60.0

        #Building a circular Gaussian beam on the 2D mesh
        x,y = mesh
        a = (np.cos(theta)**2)/(2*sigma_major**2) + (np.sin(theta)**2)/(2*sigma_minor**2)
        b = 1*np.sin(2*theta)/(4*sigma_major**2) - np.sin(2*theta)/(4*sigma_minor**2)
        c = (np.sin(theta)**2)/(2*sigma_major**2) + (np.cos(theta)**2)/(2*sigma_minor**2)
        norm_factor = 2*np.pi*sigma_major*sigma_minor
        beam_kernel = np.exp(-1*(a*(x - x0)**2 + 2*b*(x - x0)*(y - y0) + c*(y - y0)**2)) 
        beam_kernel /= norm_factor
        
        size = beam_kernel.shape[0]
        self.beam_kernel = np.zeros((4, size, size))
        self.beam_kernel[0] = beam_kernel
        self.beam_kernel[1] = -1.0*beam_kernel


    def get_beam_row(del_beta):
        row_num = np.where(self.del_beta==del_beta)[0][0]

        return self.beam_kernel[...,row_num]


    def check_normalisation(self):
        dx = self.config.beam_resolution
        dy = self.config.beam_resolution
        map_index = ["T", "Q", "U", "V"]
        for i in range(4): 
            integral = np.sum(self.beam_kernel[i])*dx*dy
            print "The integral of the", map_index[i], "beam is :", integral
            print "Percentage difference with unity :", 100*(1-integral)


    def get_mesh(self):
        offset_max = max(abs(self.config.offset_x/60.0), abs(self.config.offset_y/60.0))              #arc-mins
        fwhm_major = self.config.fwhm_major
        fwhm_minor = self.config.fwhm_minor
        size = self.config.beam_cutoff*fwhm_major + offset_max             #arc-mins
        dd = self.config.beam_resolution                                            #arc-mins
        n = int(size/dd/2)
        x = np.arange(-n, n+1)*dd
        y = -1*np.arange(-n, n+1)*dd
        self.del_beta = x
        return np.meshgrid(x,y)


    def display_beam_settings(self):
        factor = 2*np.sqrt(2*np.log(2))
        fwhm_major = self.config.fwhm_major
        fwhm_minor = self.config.fwhm_minor
        ellipticity = 100*2*(fwhm_major - fwhm_minor)/(fwhm_major + fwhm_minor)
        print "Major axis(FWHM) :", fwhm_major, "arcmins" 
        print "Minor axis(FWHM) :", fwhm_minor, "arcmins"
        print "Ellipticity :", ellipticity, "%"
        print "Center :", self.config.offset_x, config.offset_y
        print "Tilt :", self.config.beam_angle, "degrees"
        print "Pixel size :", self.config.beam_resolution, "arcmins" 
        print "Kernel width in FWHM of beam:", self.config.beam_cutoff
        print "# of pixels per FWHM (minor-axis) of beam :", fwhm_minor/self.config.beam_resolution
        print "Expected # of pixels in kernel cross-section :", int(self.config.beam_cutoff*fwhm_major/self.config.beam_resolution/2)*2 + 1 


    def plot_beam(self):
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        n = self.beam_kernel[0].shape[0]/2
        extent = np.arange(-n, n+1)*self.config.beam_resolution
        im = new_imshow(ax1, self.beam_kernel[0], x=extent, y=extent, interpolation="nearest")
        fig.colorbar(im, ax=ax1)
        im = new_imshow(ax2, self.beam_kernel[1], x=extent, y=extent, interpolation="nearest")
        fig.colorbar(im, ax=ax2)
        im = new_imshow(ax3, self.beam_kernel[2], x=extent, y=extent, interpolation="nearest")
        fig.colorbar(im, ax=ax3)
        im = new_imshow(ax4, self.beam_kernel[3], x=extent, y=extent, interpolation="nearest")
        fig.colorbar(im, ax=ax4)
        plt.show()


    def write_beam(self, out_dir=None):
        if out_dir == None:
            out_dir = self.config.out_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        np.save(os.path.join(out_dir, self.config.beam_file_name), self.beam_kernel)


if __name__=="__main__":

    config_file = sys.argv[1]
    bolo_name = sys.argv[2]
    config = importlib.import_module("simulation.timestream_simulation.config_files." + config_file).config
    bolo_config = importlib.import_module("simulation.timestream_simulation.bolo_config_files." + config.bolo_config_file).bolo_config

    bolo_beam = Beam(config, bolo_config.bolos[bolo_name])

    if config.check_normalisation:
        bolo_beam.check_normalisation()
    if config.display_beam_settings:
        bolo_beam.display_beam_settings()
    if config.plot_beam:
        bolo_beam.plot_beam()
    if config.write_beam:
        bolo_beam.write_beam()
