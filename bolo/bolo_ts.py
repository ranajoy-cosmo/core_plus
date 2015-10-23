#!/usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys, copy
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix
from pysimulators import BeamGaussian
from pysimulators.interfaces.healpy import SceneHealpixCMB
from pysimulators.interfaces.healpy import HealpixConvolutionGaussianOperator
import simulation.pointing.generate_pointing as gen_p
from simulation.beam.beam_kernel import get_beam

def Bolo:

    def __init__(self, scan_settings=None, bolo_name='0001'):
        self.bolo_parmas = np.load(os.path.join(settings.bolo_params, bolo_name, ".py"))
        if scan_settings is None:
            from local_settings import settings 
            self.scan_settings = settings
        else
            self.scan_settings = scan_settings

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Simulating the time-ordered data for a pencil beam
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def simulate_tod(settings):

        #*#*#* Generating the time-ordered pointing #*#*#*
        
        if settings.generate_pointing:
            import pointing_settings
            v = gen_p.generate_pointing(pointing_settings.settings)
        elif settings.load_pointing:
            v = np.load(settings.data_folder + "pointing_0.npy") 
        else:
            print "Not loaded any pointing. Exiting"
            sys.exit()
        
        #*#*#* Generating the time-ordered set of hit pixels from the pointing #*#*#*

        hit_pix = hp.vec2pix(settings.nside, v[...,0], v[...,1], v[...,2])
        
        #*#*#* Building the projection matrix P #*#*#*

        nsamples = hit_pix.size
        npix = hp.nside2npix(settings.nside)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32,
                               dtype_index = np.int32)
        matrix.data.index = hit_pix[..., None]
        matrix.data.value = 1
        P = ProjectionOperator(matrix, shapein = npix, shapeout = nsamples)        
        
        #*#*#* Loading the input sky map #*#*#*
        sky_map = hp.read_map(settings.input_map)
        if settings.nside is not hp.get_nside(sky_map):
            print "NSIDE of the loaded map does not match with required NSIDE"
            sys.exit()

        #*#*#* Generating the time ordered signal #*#*#*
        signal = P(sky_map)

        #*#*#* Generating the hitmap #*#*#*
        hitmap = P.T(np.ones(nsamples, dtype=np.float32))
       
        #*#*#* Generating the scan patch mask #*#*#*
        mask = hitmap > 0
        scanned_map = np.full(hitmap.size, np.nan)
        scanned_map[mask] = sky_map[mask]

        if settings.display_scanned_map:
            hp.mollzoom(hitmap)
            plt.show()
            hp.mollzoom(scanned_map)
            plt.show()
        
        if settings.write_scanned_map:
            hp.write_map(settings.output_folder + "hitmap.fits", hitmap)
            hp.write_map(settings.output_folder + "scanned_map.fits", scanned_map)

        if settings.write_signal:
            np.save(settings.data_folder + "signal", signal)
        elif settings.pipe_with_map_maker:
            return signal, P

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Simulating the time-ordered data for an extended beam
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def simulate_beam_tod(self):
        #Generating the pointing
        beam_kernel, del_beta = get_beam()
        import local_settings.pointing_params as pointing_params
        beta_0 = pointing_params.beta_0
        beta = beta_0 + del_beta

        sky_map = self.load_maps(scan_settings)
        
        #Building the projection matrix P
        nsamples = int(1000.0*pointing_params.t_flight/pointing_params.t_sampling) 
        npix = hp.nside2npix(scan_settings.nside)
        matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32,
                               dtype_index = np.int32)
        matrix.data.value = 1
        P = ProjectionOperator(matrix, shapein = npix, shapeout = nsamples)        
         
        signal = np.zeros(nsamples)
        
        for i in range(beta.size):
            v = gen_p.generate_pointing(pointing_settings.settings, beta[i])
            hit_pix = hp.vec2pix(settings.nside, v[...,0], v[...,1], v[...,2])
        
            matrix.data.index = hit_pix[..., None]
            if i is beta.size/2:
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
            np.save(settings.output_folder + "signal", signal)

        if settings.display_scanned_map:
            hp.mollzoom(hitmap)
            plt.show()
            hp.mollzoom(scanned_map)
            plt.show()
        
        if settings.write_scanned_map:
            hp.write_map(settings.output_folder + "hitmap.fits", hitmap)
            hp.write_map(settings.output_folder + "scanned_map.fits", scanned_map)

        if settings.write_signal:
            np.save(settings.output_folder + "signal", signal)
        
        if settings.pipe_with_map_maker:
            return signal, P

    def load_maps(scan_settings):
        input_map = hp.read_map(scan_settings.input_mapi, field=(0,1,2))
        return input_map


if __name__=="__main__":
    from local_settings import settings
    for bolo_name in settings.bolo_names:
        bolo = Bolo(settings, bolo_name)

