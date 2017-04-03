#!/usr/bin_env python

import numpy as np
import healpy as hp
from simulation.lib.utilities.prompter import prompt

def fill_empty_pixels(sky_map, max_iter, fail_fill_value=0, pol=True):
    if np.sum(np.isnan(sky_map)) == 0:
        return

    nside = hp.get_nside(sky_map)

    if pol:
        dim = sky_map.shape[0]
        for i in xrange(max_iter):
            empty_pix = np.where(np.isnan(sky_map[0]))[0]
            theta, phi = hp.pix2ang(nside, empty_pix)
            neighbours = hp.get_all_neighbours(nside, theta, phi).T
            for j in range(dim):
                fill_values = np.nanmean(sky_map[j][neighbours], axis=-1)
                sky_map[j][empty_pix] = fill_values
            if np.sum(np.isnan(sky_map)) == 0:
                break
    else:
        for i in xrange(max_iter):
            empty_pix = np.where(np.isnan(sky_map))[0]
            theta, phi = hp.pix2ang(nside, empty_pix)
            neighbours = hp.get_all_neighbours(nside, theta, phi).T
            fill_values = np.nanmean(sky_map[neighbours], axis=-1)
            sky_map[empty_pix] = fill_values
            if np.sum(np.isnan(sky_map)) == 0:
                break

    num_empty_pix = np.sum(np.isnan(sky_map))
    if num_empty_pix:
        prompt("{} empty pixels remaining after {} iterations. Filling empty pixels with {}\n".format(num_empty_pix, max_iter, fail_fill_value))
        sky_map[np.isnan(sky_map)] = fail_fill_value
