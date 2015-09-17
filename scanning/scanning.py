#!/usr/bin/env python

import numpy as np
import healpy as hp
from settings import settings
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix
from pysimulators import BeamGaussian
from pysimulators.interfaces.healpy import SceneHealpixCMB
from pysimulators.interfaces.healpy import HealpixConvolutionGaussianOperator
import sys
import qubic
import simulation.pointing.generate_pointing as gen_p

class Scanning:
    
    def __init__(self):
        pass

    def get_hit_pixels(self):
        v = gen_p.generate_pointing().T
        self.hit_pix = hp.vec2pix(settings.nside, v[0], v[1], v[2])

    def get_proj_op(self):
        try:
            self.hit_pix
        except AttributeError:
            self.get_hit_pixels()
        nsamples = self.hit_pix.size
        npix = hp.nside2npix(settings.nside)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32,
                           dtype_index = np.int32)
        matrix.data.index = self.hit_pix[..., None]
        matrix.data.value = 1
        self.P = ProjectionOperator(matrix, shapein = npix, shapeout = nsamples)

    def get_hitmap(self):
        try:
            self.P
        except AttributeError:
            self.get_proj_op()
        nsamples = self.hit_pix.size
        self.hitmap = (self.P).T(np.ones(nsamples, dtype = np.float32))

    def load_map(self):
        self.sky_map = hp.read_map(settings.map_file, field = (0,1,2))

    def simulate_tod(self):
        self.tod = P(sky)



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

beam_weight = get_beam_weight()

beam_weight = 1
#Making the time ordered data---------------------------------------------------
y = beam_weight*P(sky)#(sky_power) + sigma_noise*np.random.standard_normal(nsamples)
#-------------------------------------------------------------------------------

def get_map():
    #Making the map-------------------------------------------------------------
    x = (P.T * P).I * P.T * y
    #x[np.isnan(x)]=0
    #return x/beam_weight
    return x
    #---------------------------------------------------------------------------

map = get_map()

"""


