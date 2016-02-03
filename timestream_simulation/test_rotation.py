#! /usr/bin/env python 

import matplotlib.pyplot as plt
import simulation.lib.plotting.vector_plot_3d as vec
import simulation.lib.quaternion.quaternion as qt

def draw(u, q):
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    vec.draw_axes(ax, length = 1)
    vec.draw_vector(u, ax, color = 'g')
    v = qt.transform(q, u)
    if len(v.shape) == 1:
        vec.draw_vector(v, ax, color = 'r')
    else:
        for vt in v:
            vec.draw_vector(vt, ax, color = 'r')
    plt.show()
