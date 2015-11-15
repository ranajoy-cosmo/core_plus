#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
import simulation.bolo.timestream_interp as ts
from custom_settings import settings, scan_params
import os

def make_map_from_signal():

    signal = ts.do_simulation(scan_params)
    v = np.load(os.path.join(settings.output_folder, "pointing.npy"))
    hit_pix = hp.vec2pix(settings.nside_out, v[...,0], v[...,1], v[...,2])
    nsamples = hit_pix.size
    npix = hp.nside2npix(settings.nside_out)
    matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.value = 1
    matrix.data.index = hit_pix[..., None]
    P = ProjectionOperator(matrix, shapein = npix, shapeout = nsamples)
    P.matrix = matrix
    sky_map = (P.T*P).I*P.T*signal

    if settings.display_map:
        hp.mollzoom(sky_map)
        plt.show()

    hitmap = P.T(np.ones(nsamples, dtype=np.float32))

    if settings.write_map:
        hp.write_map(os.path.join(settings.output_folder, "reconstructed_map.fits"), sky_map)
        hp.write_map(os.path.join(settings.output_folder, "hitmap_out.fits"), hitmap)


if __name__=="__main__":
    make_map_from_signal()
    
