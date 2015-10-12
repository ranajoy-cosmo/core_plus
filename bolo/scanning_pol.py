#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
from pysimulators import BeamGaussian
from pysimulators.interfaces.healpy import SceneHealpixCMB
from pysimulators.interfaces.healpy import HealpixConvolutionGaussianOperator
import sys, copy
import qubic
import simulation.pointing.generate_pointing as gen_p

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for a pencil beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

class Bolo:

    def __init__(self, settings=None):
        if settings is None:
            import local_settings
            self.settings = local_settings.settings
        else:
            self.settings = settings


    def simulate_tod(self):
        #*#*#* Generating the time-ordered pointing #*#*#*
        
        if self.settings.generate_pointing:
            import pointing_settings
            pointing_settings.do_pol = self.settings.do_pol
            if self.settings.do_pol:
                v, pol_ang = gen_p.generate_pointing(pointing_settings.settings)
            else:
                v = gen_p.generate_pointing(pointing_settings.settings)
        elif self.settings.load_pointing:
            v = np.load(settings.data_folder + "pointing_0.npy") 
        else:
            print "Not loaded any pointing. Exiting"
            sys.exit()
        
        #*#*#* Generating the time-ordered set of hit pixels from the pointing #*#*#*

        hit_pix = hp.vec2pix(settings.nside, v[...,0], v[...,1], v[...,2])
        
        #*#*#* Building the projection matrix P #*#*#*

        nsamples = hit_pix.size
        npix = hp.nside2npix(self.settings.nside)
        if self.settings.do_pol:
            matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index=np.int32)
            matrix.data.index[:, 0] = hit_pix
            matrix.data.value[:, 0, 0, 0] = 0.5
            matrix.data.value[:, 0, 0, 1] = 0.5 * np.cos(2*pol_ang)
            matrix.data.value[:, 0, 0, 2] = 0.5 * np.sin(2*pol_ang)
        else:
            matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32,
                                   dtype_index = np.int32)
            matrix.data.index = hit_pix[..., None]
            matrix.data.value = 1
        P = ProjectionOperator(matrix)#, shapein = npix, shapeout = nsamples)        
        
        #Generating the hitmap
        hitmap = P.T(np.ones(nsamples, dtype=np.float32))
        if settings.do_pol:
            hitmap = hitmap[:, 0]*2
       
        #Generating the scan patch


        #*#*#* Loading the input sky map #*#*#*
        sky_map = np.array(self.load_input_map()) 
        valid = hitmap > 0
        scanned_map = np.full(sky_map.shape, np.nan)
        scanned_map[...,valid] = sky_map[...,valid]

        #Generating the time ordered signal
        if self.settings.do_pol:
            P = P.restrict(valid, inplace=True)
            pack = PackOperator(valid, broadcast='rightward')
            signal = P*pack*sky_map
        else:
            signal = P(sky_map)
        
        if self.settings.display_scanned_map:
            self.display_scanned_map(scanned_map)

        if self.settings.write_scanned_map:
            hp.write_map(self.settings.output_folder + "hitmap.fits", hitmap)
            hp.write_map(self.settings.output_folder + "scanned_map.fits", scanned_map)

        if self.settings.write_signal:
            np.save(self.settings.output_folder + "signal", signal)
        elif settings.pipe_with_map_maker:
            return signal, P

    def load_input_map(self):
        if self.settings.do_pol:
            sky_map = hp.read_map(self.settings.input_map, field = (0,1,2))
        else:
            sky_map = hp.read_map(self.settings.input_map)

        nside = hp.get_nside(sky_map)
        if nside != self.settings.nside:
            print "NSIDE of loaded map and settings do not match"
            sys.exit()
        return sky_map

    def display_scanned_map(self, scanned_map):
        if self.settings.do_pol:
            for i in (0,1,2):
                hp.mollzoom(scanned_map[i])
                plt.show()
        else:
            hp.mollzoom(scanned_map)
            plt.show()


"""

def get_conversion_factor(a):
    scene = SceneHealpixCMB(NSIDE, absolute = a)
    return scene.get_unit_conversion_operator(nu)



def get_unbeamed_absolute_sky():
    conversion = get_conversion_factor(True)
    sky_abs = conversion(T) + sky
    return sky_abs*nu*0.25

sky_power = get_unbeamed_absolute_sky()

def get_convolved():
    convolved = hp.smoothing(sky, fwhm = np.radians(beam_fwhm/60))
    #conversion = get_conversion_factor(True)
    #return conversion(T) + convolved
    return T + convolved

sky_power = get_convolved()


def get_beam_weight():
    beam = BeamGaussian(np.radians(beam_fwhm*20/60.0))
    beam_healpix = beam.healpix(NSIDE)
    return np.sum(beam_healpix)

"""

if __name__=="__main__":
    from local_settings import settings
    if settings.do_beam_kernel:
        simulate_beam_tod(settings)
    else:
        bolo = Bolo(settings)
        bolo.simulate_tod()
