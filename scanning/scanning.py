#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from settings import settings
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix
from pysimulators import BeamGaussian
from pysimulators.interfaces.healpy import SceneHealpixCMB
from pysimulators.interfaces.healpy import HealpixConvolutionGaussianOperator
import sys
import qubic
import simulation.pointing.generate_pointing as gen_p
from beam_kernel import beam_kernel

def simulate_tod():
    #Generating the pointing
    if settings.generate_pointing:
        gen_p.generate_pointing()
    v = np.load(settings.data_folder + "pointing_0.npy") 
    hit_pix = hp.vec2pix(settings.nside, v[...,0], v[...,1], v[...,2])
    
    #Building the projection matrix P
    nsamples = hit_pix.size
    npix = hp.nside2npix(settings.nside)
    matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32,
                           dtype_index = np.int32)
    matrix.data.index = hit_pix[..., None]
    matrix.data.value = 1
    P = ProjectionOperator(matrix, shapein = npix, shapeout = nsamples)        
    
    sky_map = hp.read_map(settings.input_map_path)

    #Generating the time ordered signal
    signal = P(sky_map)
    signal = np.convolve(signal, beam_kernel, mode = 'same')

    #Generating the hitmap
    hitmap = P.T(np.ones(nsamples, dtype=np.float32))
   
    #Generating the scan patch
    mask = np.full(hitmap.size, False, dtype = bool)
    mask[hitmap > 0] = True
    scanned_map = np.full(hitmap.size, np.nan)
    scanned_map[mask] = sky_map[mask]

    if settings.display_scanned_map:
        hp.mollzoom(hitmap)
        plt.show()
        hp.mollzoom(scanned_map)
        plt.show()
    
    if settings.write_scanned_map:
        hp.write_map(settings.data_folder + "hitmap.fits", hitmap)
        hp.write_map(settings.data_folder + "scanned_map.fits", scanned_map)

    if settings.write_signal:
        np.save(settings.data_folder + "signal", signal)
    elif settings.pipe_with_map_maker:
        return signal, P

    

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


