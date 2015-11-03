#!/usr/bin/env python

import numpy as np
import healpy as hp
import rotation_pol as r
from params import *
from pysimulators import ProjectionOperator, BeamGaussian
from pysimulators.sparse import FSRBlockMatrix
from pysimulators.interfaces.healpy import SceneHealpixCMB
from pyoperators import DiagonalOperator, PackOperator, pcg, Rotation3dOperator

import sys
import qubic


def do_scan():
    #Printing the parameters for the simulation----------------------------------------------------------------------
    print "theta_cross = ", theta_cross, " arcmin"
    print "theta_co = ", theta_co, " arcmin"
    print "T_orbit = ", t_orbit/3600.0/24.0, " days"
    print "T_prec = ", t_prec, "seconds"
    print "T_spin = ", t_spin, "seconds"
    print "T_sampling = ", dt*1000, "milli-seconds"
    print "alpha = ", np.degrees(alpha), "degrees"
    print "beta = ", np.degrees(beta), "degrees"
    print "Size of pix array ", 8.0*t_orbit/(1000000*dt*split), "MB"
    print "NSIDE = ", NSIDE
    print "Size of hitmap array ", 8.0*npix/1000000, "MB"
    #----------------------------------------------------------------------------------------------------------------

    o = raw_input("\nDo you want to continue (y/n)? ")
    if o is 'y':pass
    else:sys.exit()

    
    time = dt*np.arange(int(t_orbit/dt/split))
    
    v, pol_dir = r.rotate(time, w_orbit, w_prec, alpha, w_spin, v_init, eta)
    v=v.T
    
    return hp.vec2pix(NSIDE,v[0],v[1],v[2]), pol_dir
 
pix, pol_dir = do_scan()

def get_proj_op():
    #Creating the Projection operator P-----------------------------------------------------------------------------
    matrix = FSRBlockMatrix((nsamples, npix*3), (1, 3), ncolmax=1, dtype=np.float32, dtype_index=np.int32)
    matrix.data.index[:, 0] = pix
    matrix.data.value[:, 0, 0, 0] = 0.5
    matrix.data.value[:, 0, 0, 1] = 0.5 * np.cos(2*pol_dir)
    matrix.data.value[:, 0, 0, 2] = 0.5 * np.sin(2*pol_dir)
    return ProjectionOperator(matrix)
    #---------------------------------------------------------------------------------------------------------------

P = get_proj_op()


#Making the hitmap-----------------------------------------------------------------------------------------------
hitmap = P.T(np.ones(nsamples, dtype=np.float32))[:, 0] * 2
mask = hitmap > 0

P = P.restrict(mask, inplace=True)
pack = PackOperator(mask, broadcast='rightward')
#----------------------------------------------------------------------------------------------------------------

def get_conversion_factor(a):
    scene = SceneHealpixCMB(NSIDE, absolute = a)
    return scene.get_unit_conversion_operator(nu)



def get_sky():
    #Simulating the CMB Sky------------------------------------------------------------------------------------------
    spectra = qubic.read_spectra(0)
    sky_del = hp.synfast(spectra, NSIDE)
    #conversion = get_conversion_factor(False)
    #return conversion(sky_del)
    return sky_del
    #----------------------------------------------------------------------------------------------------------------

sky = get_sky()

def get_convolved():
    convolved = hp.smoothing(sky, fwhm = np.radians(beam_fwhm/60.0), pol = True)
    #conversion = get_conversion_factor(True)
    #return conversion(T) + convolved
    return np.array(convolved)

sky_power = get_convolved().T

def get_beam_weight():
    beam = BeamGaussian(np.radians(beam_fwhm/60.0))
    beam_healpix = beam.healpix(NSIDE)
    return np.sum(beam_healpix)

beam_weight = get_beam_weight()

#Making the time ordered data------------------------------------------------------------------------------------
y = beam_weight*P*pack*(sky_power) #+ sigma_noise*np.random.standard_normal(nsamples)
#----------------------------------------------------------------------------------------------------------------

def get_map():
    #Making the map--------------------------------------------------------------------------------------------------
    A = P.T * P
    b = P.T * y
    M = DiagonalOperator(1 / hitmap[mask], broadcast='rightward')
    # solve for Ax = b
    solution = pcg(A, b, M=M, disp=True, tol=3*1e-4, maxiter=2000)
    x = pack.T * solution['x']
    x[hitmap == 0] = np.nan
    return x/beam_weight
    #----------------------------------------------------------------------------------------------------------------

map = get_map().T
