#! /usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.ndimage import zoom, interpolation
import sys
import os
import importlib
from simulation.lib.plotting.my_imshow import new_imshow
from simulation.lib.utilities.generic_class import Generic

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
        beam_kernel = np.array(hp.mrdfits(self.config.input_beam_file))
        mark_orig_dim = np.sqrt(beam_kernel.shape[1])                             #pixels
        beam_kernel = beam_kernel.reshape((4, mark_orig_dim, mark_orig_dim))
        mark_fwhm_major = 7.68                   #arc-mins
        mark_fwhm_minor = 7.68                   #arc-mins
        mark_resolution = 0.39523370660946627            #arc-mins

        new_dim = int(mark_orig_dim * mark_resolution / self.config.scan_resolution)
        if new_dim%2 == 0:
            new_dim += 1
        #print "old resolution :", mark_resolution
        #print "new_resolution :", self.config.scan_resolution
        old_extent = mark_orig_dim * mark_resolution
        new_extent = new_dim * self.config.scan_resolution
        #print "old extent :", mark_resolution * 181
        #print "new_extent :", new_extent
        self.config.fwhm_major = mark_fwhm_major * new_extent / old_extent
        self.config.fwhm_minor = mark_fwhm_minor * new_extent / old_extent
        #print "old fwhm :", mark_fwhm_major
        #print "new fwhm :", self.config.fwhm_major

        self.beam_kernel = np.empty((4, new_dim, new_dim))
        for i in range(4):
            self.beam_kernel[i] = zoom(beam_kernel[i], float(new_dim)/float(mark_orig_dim))
            self.beam_kernel[i] = interpolation.rotate(self.beam_kernel[i], angle=self.config.beam_angle, reshape=False)

        #print "beam_cutoff required :", self.config.beam_cutoff
        beam_extension_required = self.config.fwhm_major * self.config.beam_cutoff          #arc-min
        #print "extention required :", beam_extension_required
        num_pix = int(new_dim * beam_extension_required / new_extent)
        if num_pix%2 == 0:
            num_pix -= 1
        #print "num pix :", num_pix
        #print "extension achieved :", num_pix * self.config.scan_resolution

        self.config.beam_cutoff = num_pix * self.config.scan_resolution / self.config.fwhm_major
        #print "cutoff achieved :", self.config.beam_cutoff

        start = new_dim/2 - num_pix/2
        stop = new_dim/2 + num_pix/2 + 1
        self.beam_kernel = self.beam_kernel[..., start:stop, start:stop]
        #print "kernel shape :", self.beam_kernel.shape
        self.del_beta = self.config.scan_resolution * np.arange(-num_pix/2 + 1, num_pix/2 + 1)


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


    def get_beam_row(self, del_beta):
        row_num = np.where(self.del_beta==del_beta)[0][0]

        return self.beam_kernel[...,row_num]


    def check_normalisation(self):
        dx = self.config.scan_resolution
        dy = self.config.scan_resolution
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
        dd = self.config.scan_resolution                                            #arc-mins
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
        print "Center :", self.config.offset_x, self.config.offset_y
        print "Tilt :", self.config.beam_angle, "degrees"
        print "Beam pixel size :", self.config.scan_resolution, "arcmins" 
        print "Kernel width in FWHM of beam:", self.config.beam_cutoff
        print "# of pixels per FWHM (minor-axis) of beam :", fwhm_minor/self.config.scan_resolution
        sys.stdout.flush()


    """
    def plot_beam(self):
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        n = self.beam_kernel[0].shape[0]/2
        extent = np.arange(-n, n+1)*self.config.scan_resolution
        im = new_imshow(ax1, self.beam_kernel[0], x=extent, y=extent, interpolation="nearest")
        ax1.set_title('T')
        fig.colorbar(im, ax=ax1)
        im = new_imshow(ax2, self.beam_kernel[1], x=extent, y=extent, interpolation="nearest")
        ax2.set_title('Q')
        fig.colorbar(im, ax=ax2)
        im = new_imshow(ax3, self.beam_kernel[2], x=extent, y=extent, interpolation="nearest")
        ax3.set_title('U')
        fig.colorbar(im, ax=ax3)
        im = new_imshow(ax4, self.beam_kernel[3], x=extent, y=extent, interpolation="nearest")
        ax4.set_title('V')
        fig.colorbar(im, ax=ax4)
        fig.suptitle("Rescaled Plack 217_5a, FWHM : 7.68', Resol : 0.96', Extent : 3.85*FWHM")
        plt.show()
    """
    def plot_beam(self):
        fig, ax= plt.subplots()
        n = self.beam_kernel[0].shape[0]/2
        extent = np.arange(-n, n+1)*self.config.scan_resolution
        im = new_imshow(ax, 10*self.beam_kernel[0], x=extent, y=extent, interpolation="nearest", cmap='gray')
        plt.show()


    def write_beam(self, out_dir=None):
        if out_dir == None:
            out_dir = os.path.join(os.getcwd(), "beam_maps") 
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        np.save(os.path.join(out_dir, self.config.beam_file_name), self.beam_kernel)


if __name__=="__main__":

    config_file = sys.argv[1]
    bolo_name = sys.argv[2]
    config = importlib.import_module("simulation.timestream_simulation.config_files." + config_file).config
    config.scan_resolution = 1.5 #0.959620009804 
    bolo_config = importlib.import_module("simulation.timestream_simulation.bolo_config_files." + config.bolo_config_file).bolo_config

    bolo_beam = Beam(config, bolo_config.bolos[bolo_name])

    #if config.check_normalisation:
    #    bolo_beam.check_normalisation()
    #if config.display_beam_settings:
    #    bolo_beam.display_beam_settings()
    bolo_beam.plot_beam()
    #if config.write_beam:
    #    bolo_beam.write_beam()
