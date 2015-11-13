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

    def __init__(self, settings=None, pointing_params=None, beam_params=None, bolo_name='0001'):
        if pointing_params is None:
            from custom_settings import pointing_params 
            self.pointing_params = pointing_params  
        else:
            self.pointing_params = pointing_params

        if beam_params is None:
            from custom_settings import beam_params 
            self.beam_params = beam_params  
        else:
            self.beam_params = beam_params

        if settings is None:
            from custom_settings import settings 
            self.settings = settings
        else:
            self.settings = settings


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Simulating the time-ordered data for any beam
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def simulate_timestream(self):
        
        #Getting the beam profile and the del_beta
        beam_kernel, del_beta = get_beam()

        sky_map = hp.read_map(self.settings.input_map)

        nsamples = int(1000.0*self.pointing_params.t_flight/self.pointing_params.t_sampling)*self.settings.oversampling_rate 
        signal = np.zeros(nsamples)
        
        for i in range(del_beta.size):
            v = gen_p.generate_pointing(self.pointing_params, np.deg2rad(del_beta[i]/60.0))
            theta, phi = hp.vec2ang(v)
        
            if i is del_beta.size/2: 
                v_central = v[::self.settings.oversampling_rate]
                
            #Generating the time ordered signal
            signal += np.convolve(hp.get_interp_val(sky_map, theta, phi), beam_kernel.T[i], mode = 'same')

        beam_sum = np.sum(beam_kernel)
        signal/=beam_sum
        signal = signal[::self.settings.oversampling_rate]

        hitmap = self.get_hitmap(v_central)

        scanned_map = self.get_scanned_map(sky_map, hitmap)

        if self.settings.write_signal:
            np.save(os.path.join(self.settings.output_folder, "signal"), signal)

        if self.settings.display_scanned_map:
            hp.mollzoom(hitmap)
            plt.show()
            hp.mollzoom(scanned_map)
            plt.show()
        
        if self.settings.write_scanned_map:
            hp.write_map(os.path.join(self.settings.output_folder, "hitmap.fits"), hitmap)
            hp.write_map(os.path.join(self.settings.output_folder, "scanned_map.fits"), scanned_map)

        if self.settings.write_signal:
            np.save(os.path.join(self.settings.output_folder, "signal"), signal)
        
        if self.settings.pipe_with_map_maker:
            return signal

    def get_hitmap(self, v):
        nsamples = v.size
        npix = hp.nside2npix(self.settings.nside_in)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
        matrix.data.value = 1
        P = ProjectionOperator(matrix, shapein=npix, shapeout=nsamples)        
        hit_pix = hp.vec2pix(self.settings.nside_in, v[...,0], v[...,1], v[...,2])
        matrix.data.index = hit_pix[..., None]
        P.matrix = matrix
        hitmap = P.T(np.ones(nsamples, dtype=np.float32))
        return hitmap

    def get_scanned_map(self, sky_map, hitmap=None):
        if hitmap is None:
            print "Need hitmap first"
            sys.exit()
        mask = np.full(hitmap.size, False, dtype = bool)
        mask[hitmap > 0] = True
        scanned_map = np.full(hitmap.size, np.nan)
        scanned_map[mask] = sky_map[mask]
        return scanned_map

    def load_maps(self):
        input_map = hp.read_map(self.settings.input_map, field=(0,1,2))
        return input_map

def do_simulation(settings=None):
    if settings==None:
        from custom_settings import settings
    bolo = Bolo(settings=settings)
    if settings.pipe_with_map_maker:
        return bolo.simulate_timestream()

if __name__=="__main__":
    from custom_settings import settings
    for bolo_name in settings.bolo_names:
        bolo = Bolo()
        bolo.simulate_timestream()
