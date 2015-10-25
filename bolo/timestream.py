#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys, copy, os
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix
from pysimulators import BeamGaussian
from pysimulators.interfaces.healpy import SceneHealpixCMB
from pysimulators.interfaces.healpy import HealpixConvolutionGaussianOperator
import simulation.pointing.generate_pointing as gen_p
from simulation.beam.beam_kernel import get_beam

class Bolo:

    def __init__(self, settings=None, pointing_params=None, bolo_name='0001'):
        if pointing_params is None:
            from local_settings import pointing_params 
            self.pointing_params = pointing_params  
        else:
            self.pointing_params = pointing_params

        if settings is None:
            from local_settings import settings 
            self.settings = settings
        else:
            self.settings = settings

        #self.bolo_params = np.load(os.path.join(settings.bolo_param_folder, bolo_name) + '.py')

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def simulate_timestream(self):
        
        #Getting the beam profile and the del_beta
        beam_kernel, del_beta = get_beam()

        sky_map = hp.read_map(settings.input_map)

        #Building the projection matrix P
        nsamples = int(1000.0*self.pointing_params.t_flight/self.pointing_params.t_sampling) 
        npix = hp.nside2npix(self.settings.nside)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32,
                               dtype_index = np.int32)
        matrix.data.value = 1
        P = ProjectionOperator(matrix, shapein = npix, shapeout = nsamples)        
         
        signal = np.zeros(nsamples)
        
        for i in range(del_beta.size):
            v = gen_p.generate_pointing(self.pointing_params, del_beta[i])
            hit_pix = hp.vec2pix(self.settings.nside, v[...,0], v[...,1], v[...,2])
        
            matrix.data.index = hit_pix[..., None]
            if i is del_beta.size/2:
                matrix_central = copy.deepcopy(matrix)
            P.matrix = matrix
            #Generating the time ordered signal
            signal += np.convolve(P(sky_map), beam_kernel.T[i], mode = 'same')

        beam_sum = np.sum(beam_kernel)
        signal/=beam_sum

        P.matrix = matrix_central
        
        #Generating the hitmap
        hitmap = P.T(np.ones(nsamples, dtype=np.float32))
       
        #Generating the scan patch
        mask = np.full(hitmap.size, False, dtype = bool)
        mask[hitmap > 0] = True
        scanned_map = np.full(hitmap.size, np.nan)
        scanned_map[mask] = sky_map[mask]

        if settings.write_signal:
            np.save(self.settings.output_folder + "signal", signal)

        if self.settings.display_scanned_map:
            hp.mollzoom(hitmap)
            plt.show()
            hp.mollzoom(scanned_map)
            plt.show()
        
        if self.settings.write_scanned_map:
            hp.write_map(self.settings.output_folder + "hitmap.fits", hitmap)
            hp.write_map(self.settings.output_folder + "scanned_map.fits", scanned_map)

        if self.settings.write_signal:
            np.save(self.settings.output_folder + "signal", signal)
        
        if self.settings.pipe_with_map_maker:
            return signal, P

    def load_maps(self):
        input_map = hp.read_map(self.settings.input_mapi, field=(0,1,2))
        return input_map


if __name__=="__main__":
    from local_settings import settings
    for bolo_name in settings.bolo_names:
        bolo = Bolo()
        bolo.simulate_timestream()
