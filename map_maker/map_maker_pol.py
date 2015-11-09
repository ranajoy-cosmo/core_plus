#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from pyoperators import DiagonalOperator, PackOperator, pcg
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
import simulation.bolo.timestream_pol as ts
from custom_settings import settings, scan_params
import os

def make_map_from_signal():

    signal = ts.do_simulation(scan_params)
    v = np.load(os.path.join(settings.output_folder, "pointing.npy"))
    pol_ang = np.load(os.path.join(settings.output_folder, "pol_angle.npy"))

    hit_pix = hp.vec2pix(settings.nside_out, v[...,0], v[...,1], v[...,2])
    nsamples = hit_pix.size
    npix = hp.nside2npix(settings.nside_out)
    
    matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.value = 1
    matrix.data.index[:, 0] = hit_pix
    matrix.data.value[:, 0, 0, 0] = 0.5
    matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol_ang)
    matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol_ang)
    
    P = ProjectionOperator(matrix)#, shapein = npix, shapeout = nsamples)
    
    hitmap = P.T(np.ones(nsamples, dtype=np.float32))[:, 0]*2
    mask = hitmap>0

    P = P.restrict_new(mask, inplace=True)
    pack = PackOperator(mask, broadcast='rightward')

    A = P.T*P
    b = P.T*signal
    M = DiagonalOperator(1/hitmap[mask], broadcast='rightward')

    solution = pcg(A, b, M=M, disp=True, tol=1e-4, maxiter=2000)
    x = pack.T*solution['x']
    x[hitmap == 0] = np.nan
    sky_map = x.T

    if settings.write_map:
        hp.write_map(os.path.join(settings.output_folder, "reconstructed_map.fits"), sky_map)
        hp.write_map(os.path.join(settings.output_folder, "hitmap_out.fits"), hitmap)

    if settings.display_map:
        hp.mollzoom(sky_map[0])
        plt.show()



if __name__=="__main__":
    make_map_from_signal()
    
