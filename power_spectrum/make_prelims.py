#!/usr/bin/env python 

import numpy as np
import pyfits as pf
import healpy as hp
import sys
import os


def make_beam():
    ell = np.arange(0,lmax + 500)
    sigma = np.radians(fwhm/60.0)/2.35482
    beam = np.exp(-0.5*ell*(ell + 1)*sigma**2)
    pf.writeto(os.path.join(map_dir, "beam.fits"), beam)

def make_bins(width=10):
    bins = np.arange(0, lmax, width)
    bins = bins.astype(float)
    pf.writeto(os.path.join(map_dir, "bins.fits"), bins)

def make_binary_mask():
    sky = hp.read_map(map_file)
    binary_hitmap = np.logical_not(np.isnan(sky))
    hp.write_map(os.path.join(map_dir, "binary_mask.fits"), binary_hitmap)
    hp.write_map(os.path.join(map_dir, "weight_map.fits"), binary_hitmap)


if __name__=="__main__":
    map_dir, map_file, lmax, fwhm = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
    make_bins(10)
    make_beam()
    make_binary_mask()
