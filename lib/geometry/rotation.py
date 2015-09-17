#!/usr/bin/env python

import numpy as np
from numpy import cos,sin

class Rotation:

    def __init__(self):
        self.R = np.empty((3,3))

    """

    def Rx(theta):
        a=np.cos(theta)
        b=np.sin(theta)
        return np.array([[1.0,0,0],[0,a,b],[0,-b,a]])

    def Ry(theta):
        a=np.cos(theta)
        b=np.sin(theta)
        return np.array([[a,0,-b],[0,1.0,0],[b,0,a]])

    def Rz(theta):
        a=np.cos(theta)
        b=np.sin(theta)
        return np.array([[a,b,0],[-b,a,0],[0,0,1]])

    def rotate(phi,alpha,psi,beta):
        R1=Ry(beta)
        R2=Rx(-phi)
        R3=Ry(alpha)
        R4=Rx(-psi)

        R=np.dot(R1,np.dot(R2,np.dot(R3,R4)))

        return R

    """

    def rotate(time,w_prec,alpha,w_spin,v,eta):
        cp=cos(w_prec*time)
        sp=sin(w_prec*time)
        ct=cos(np.full(time.size, alpha))
        st=sin(np.full(time.size, alpha))
        cs=cos(w_spin*time)
        ss=sin(w_spin*time)

        R=np.array([ [ ct , -sp*st , st*cp ] ,
                    [ -st*ss , cp*cs - sp*ct*ss , sp*cs + cp*ct*ss ] ,
                     [ -st*cs , -cp*ss - sp*ct*cs , -sp*ss + cp*ct*cs ] ])

        #pol_dir=(eta+w_prec*time+w_spin*time)%np.pi

        return R, np.dot(R.T,v)#,pol_dir


    

