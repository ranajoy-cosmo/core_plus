#! /usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import pyoperators as po
import simulation.lib.plotting.vector_plot_3d as vec

"""
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
"""

def generate_pointing(settings = None):
    if settings == None:
        from settings import settings
    u_init = np.array([np.cos(settings.beta), 0.0, np.sin(settings.beta)])
    n_slices = int(1000*settings.t_flight/settings.t_sampling)
    print "No. of time steps : ", n_slices
    print "~Memory usage : ", 10*n_slices*8.0/1024/1024/1024, " GB" 
    t_steps = 0.001*settings.t_sampling*np.arange(n_slices)
    w_prec = 2*np.pi/settings.t_prec
    w_spin = 2*np.pi/settings.t_spin
    R = po.Rotation3dOperator("XY'X''", w_prec*t_steps, -1.0*np.full(n_slices, settings.alpha), w_spin*t_steps)
    v = R*u_init
    return v

if __name__ == "__main__":
    from settings import settings
    generate_pointing(settings)

