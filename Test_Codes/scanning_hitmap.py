#!/usr/bin/env python

import numpy as np
import healpy as hp
import rotation_pol as r
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix
import sys


#Mathematical constants----------------
pi=np.pi
#--------------------------------------

#Pre defined values on which the rest will depend----------------------------------------------------------------
#theta_cross=raw_input("Cross scanning rate : ")                  #arcminutes
#theta_co=raw_input("Co-scanning rate : ")                    #arcminutes
theta_cross = 5                            #arcminutes
theta_co = 5                               #arcminutes
t_orbit = 366*24*60*60.0         #seconds
t_prec = 4*24*60*60.0              #seconds
alpha = np.radians(50.0)              #radians
beta = np.radians(60.0)               #radians
splits = 200
NSIDE = 512
eta = 0
#----------------------------------------------------------------------------------------------------------------

#Derived quantities----------------------------------------------------------------------------------------------
t_spin = (1 / (15 * np.sin(alpha) ) ) * theta_cross * t_prec / (24 * 60)
dt = (1 / (5400 * np.sin(alpha) * np.sin(beta) ) ) * theta_co * theta_cross * t_prec / (24 * 60 * 60)
#t_spin = 60.0
#dt = 0.1
w_orbit = 2*pi/t_orbit
w_prec = 2*pi/t_prec
w_spin = 2*pi/t_spin
v_init = np.array([np.cos(beta), 0, np.sin(beta)])
npix = hp.nside2npix(NSIDE)
nsamples = int(t_orbit/(dt*splits))
#----------------------------------------------------------------------------------------------------------------


#Printing the parameters for the simulation----------------------------------------------------------------------
print "theta_cross = ", theta_cross, " arcmin"
print "theta_co = ", theta_co, " arcmin"
print "NSIDE = ", NSIDE
print "NSIDE resolution = ", hp.nside2resol(NSIDE, arcmin = True)
print "T_orbit = ", t_orbit/24/60/60.0, " days"
print "T_prec = ", t_prec/60/60.0, " hours"
print "T_spin = ", t_spin, " seconds"
print "T_sampling = ", dt, " seconds"
print "Sampling Frequency = ", 1/dt, " Hz"
print "alpha = ", alpha
print "beta = ", beta
print "Splits = ", splits
print "Size of time array split", 8.0 * t_orbit / (1000000 * splits * dt), "MB"
print "Size of hitmap array ", 8.0*npix/1000000, "MB"
#----------------------------------------------------------------------------------------------------------------

o=raw_input("\nDo you want to continue (y/n)? ")
if o is 'y':pass
else:sys.exit()

#Initialising the maps-------------------------------------------------------------------------------------------
hitmap=np.zeros(npix)
#m=np.zeros(npix)
#----------------------------------------------------------------------------------------------------------------

#Loading the maps------------------------------------------------------------------------------------------------
#planck = hp.ud_grade(hp.read_map('planck_70.fits',(0,1,2)),NSIDE)
#----------------------------------------------------------------------------------------------------------------

p1=[]
p2=[]
p3=[]
p4=[]
p5=[]
p6=[]

for i in xrange(100):
    print 'Running split ', i+1, ' of ', splits
    time = i*(t_orbit/splits) + dt*np.arange(int(t_orbit/(dt*splits)))

    v, pol_dir = r.rotate(time, w_orbit, w_prec, alpha, w_spin, v_init, eta)

    v = v.T
    
    pix = hp.vec2pix(NSIDE, v[0], v[1], v[2])

    p1.append(pol_dir[pix == hp.ang2pix(NSIDE,0,0)])
    p2.append(pol_dir[pix == hp.ang2pix(NSIDE,np.radians(20.0),0)])
    p3.append(pol_dir[pix == hp.ang2pix(NSIDE,np.pi/4,0)])
    p4.append(pol_dir[pix == hp.ang2pix(NSIDE,np.radians(75.0),0)])
    p5.append(pol_dir[pix == hp.ang2pix(NSIDE,np.pi/2,0)])
    p6.append(pol_dir[pix == hp.ang2pix(NSIDE,np.radians(80),0)])
    """   
    #Creating the Projection operator P---------------------------------------------------------------------------
    matrix = FSRMatrix((nsamples, npix), ncolmax=1, dtype=np.float32,
                       dtype_index=np.int32)
    matrix.data.index = pix[..., None]
    matrix.data.value = 1
    P = ProjectionOperator(matrix, shapein=npix, shapeout=nsamples)
    #--------------------------------------------------------------------------------------------------------------

    hitmap += P.T(np.ones(nsamples, dtype=np.float32))
    """
    #Making the map------------------------------------------------------------------------------------------------
    #y = P(planck[0])
    #x = (P.T * P).I * P.T * y
    #x[np.isnan(x)]=0
    #m+=x
    #--------------------------------------------------------------------------------------------------------------


