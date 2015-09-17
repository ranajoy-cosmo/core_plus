#!/usr/bin/env python

import numpy as np
from numpy import cos,sin

"""

def Rx(theta):
    a=np.cos(theta)
    b=np.sin(theta)
    one = np.ones(theta.size)
    zero = np.zeros(theta.size)
    return np.array([[one,zero,zero],[zero,a,b],[zero,-b,a]]).swapaxes(0,2)

def Ry(theta):
    a=np.cos(theta)
    b=np.sin(theta)
    one = np.ones(theta.size)
    zero = np.zeros(theta.size)
    return np.array([[a,zero,-b],[zero,one,zero],[b,zero,a]]).swapaxes(0,2)

def Rz(theta):
    a=np.cos(theta)
    b=np.sin(theta)
    one = np.ones(theta.size)
    zero = np.zeros(theta.size)
    return np.array([[a,b,zero],[-b,a,zero],[zero,zero,one]]).swapaxes(0,2)


def rotate(phi,alpha,psi,beta):
    R1=Ry(beta)
    R2=Rx(-phi)
    R3=Ry(alpha)
    R4=Rx(-psi)

    R=np.dot(R1,np.dot(R2,np.dot(R3,R4)))

    return R

"""

def rotate(time, w_year, w_prec, alpha, w_spin, v, eta = 0):

    ce = cos(w_year*time)
    se = sin(w_year*time)
    cp = cos(w_prec*time)
    sp = sin(w_prec*time)
    ct = cos(np.full(time.size, alpha))
    st = sin(np.full(time.size, alpha))
    cs = cos(w_spin*time)
    ss = sin(w_spin*time)

    #R=np.array([ [ ct , -sp*st , st*cp ] ,
    #             [ -st*ss , cp*cs - sp*ct*ss , sp*cs + cp*ct*ss ] ,
    #             [ -st*cs , -cp*ss - sp*ct*cs , -sp*ss + cp*ct*cs ] ])

    R = np.array([ [ce*ct + se*sp*st, se*ct - ce*sp*st, st*cp],
                   [-ce*st*ss - se*cp*cs + se*sp*ct*ss, -se*st*ss + ce*cp*cs - ce*sp*ct*ss, sp*cs + cp*ct*ss],
                   [-ce*st*cs + se*cp*ss + se*sp*ct*cs, -se*st*cs - ce*cp*ss - ce*sp*ct*cs, -sp*ss + cp*ct*cs] ])

    #pol_dir=(eta+w_prec*time+w_spin*time)%np.pi

    return np.dot(R.T,v)#, pol_dir


    

