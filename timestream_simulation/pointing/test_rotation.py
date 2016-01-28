#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyoperators as po
import simulation.lib.plotting.vector_plot_3d as vec
import sys

def draw(u, R):
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    vec.draw_axes(ax, length = 1)
    v = R*u
    vec.draw_vector(u, ax, color = 'g')
    if len(v.shape) == 1:
        vec.draw_vector(v, ax, color = 'r')
    else:
        for vt in v:
            vec.draw_vector(vt, ax, color = 'r')
    plt.show()
