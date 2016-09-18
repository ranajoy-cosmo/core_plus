#!/usr/bin/env python

import numpy as np
import healpy as hp
from mpi4py import MPI
from pyoperators import DiagonalOperator, PackOperator, pcg, MPIDistributionIdentityOperator
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix, FSRBlockMatrix
import os
import sys


def map_maker_function(signal_diff, v, pol_ang, grad_T):
    hit_pix = hp.vec2pix(map_making_params.nside_out, v[...,0], v[...,1], v[...,2])
    nsamples = hit_pix.size
    npix = hp.nside2npix(config.nside_out)

    matrix_A = FSRBlockMatrix((nsamples, npix*2), (1, 2), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.index[:, 0] = hit_pix
    matrix.data.value[:, 0, 0, 0] = np.cos(2*pol_ang)
    matrix.data.value[:, 0, 0, 1] = np.sin(2*pol_ang)

    matrix_grad_T = 
