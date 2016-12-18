import numpy as np
import healpy as hp

def passive_rotation(Q_in, U_in, theta, degrees=False):
    """
    Transformation of (Q,U) Stokes parameters under a rotation of the coordinate system in the counter-clockwise direction.
    Q_in and U_in are the Stokes parameters in the old system and theta are the angles by which the coordinate system rotates.
    """

    if degrees:
        theta = np.radians(theta)

    Q_out = Q_in*np.cos(2*theta) + U_in*np.sin(2*theta)
    U_out = -Q_in*np.sin(2*theta) + U_in*np.cos(2*theta)

    return Q_out, U_out
