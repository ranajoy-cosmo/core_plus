#!/usr/bin/env python

import numpy as np
import healpy as hp
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRBlockMatrix

nside = 512
npix = 12*nside**2

nsamples = 12*60*60*200 

sky = np.empty((3, npix))
sky[0] = 1.0*np.arange(npix)
sky[1] = 2.0*np.arange(npix)
sky[2] = 3.0*np.arange(npix)

hitpix = np.random.random_integers(0, npix-1, nsamples) 
pol_ang = 2*np.pi*np.random.random(nsamples)

def get_projection_matrix():

    matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index = np.int32)
    matrix.data.index[:, 0] = hitpix
    matrix.data.value[:, 0, 0, 0] = 0.5
    matrix.data.value[:, 0, 0, 1] = 0.5*np.cos(2*pol_ang) 
    matrix.data.value[:, 0, 0, 2] = 0.5*np.sin(2*pol_ang)

    return ProjectionOperator(matrix)

def get_hitmap(P):
    
    hitmap = P.T(np.ones(nsamples, dtype=np.float32))[:, 0]*2
    return hitmap

def get_signal(P, sky):

    d = P(sky.T)
    return d

def get_b_matrix(P, d):

    b = P.T(d)
    return b

def get_cov_matrix(P):

    matrix = np.squeeze(P.matrix.data.value)
    matrix_2 = np.einsum('ij...,i...->ij...', matrix, matrix)
    cov_matrix_inv = np.zeros((npix, 3, 3))
    for i in xrange(nsamples):
        cov_matrix_inv[hitpix[i]] += matrix_2[i]
    
    det = np.linalg.det(cov_matrix_inv)
    bad_pixel = (det==0.0)
    cov_matrix_inv[bad_pixel] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    return np.linalg.inv(cov_matrix_inv), bad_pixel

def get_map(cov_matrix, b):

    sky_rec = np.sum(cov_matrix*b[..., None], axis=1)
    return sky_rec.T

proj_matrix = get_projection_matrix()

signal = get_signal(proj_matrix, sky)

hitmap = get_hitmap(proj_matrix)

b = get_b_matrix(proj_matrix, signal)

cov_matrix, bad_pix = get_cov_matrix(proj_matrix)

sky_rec = get_map(cov_matrix, b)
sky_rec[..., bad_pix] = np.nan
