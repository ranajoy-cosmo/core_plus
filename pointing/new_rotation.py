#!/usr/bin/env python 

import numpy as np
import healpy as hp

pi = np.pi

NSIDE = 64

theta_in = pi/4
phi_in = 0
alpha = pi/2

PHI = np.linspace(0, pi, 180)

theta = alpha*np.cos(PHI)
phi = alpha*np.sin(PHI)

hitpix = hp.ang2pix(NSIDE, pi/2 - theta, phi)
