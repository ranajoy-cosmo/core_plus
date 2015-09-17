#!/usr/bin/env python

import math
from vector import Vector
from quaternion import Quaternion

class Rotation:
    def __init__(self,q=None,theta=None,d=None,mode=1):
        if q is None:
            if mode is 1:
                self.q=Quaternion(q_v=d,norm=1,theta=theta/2)
            else:
                self.q=Quaternion(q_v=d,norm=1,theta=-theta/2)
        else:
            self.q=q

    def rotate(self,v):
        res=self.q*Quaternion(0,v)*self.q.conjugate()
        return res.q

    def compose(*args):
        res=Quaternion(1.0,Vector(0,0,0))
        for i in args[::-1]:
            res=res*i.q
        return Rotation(res)
    
    def inverse(self):
        return Rotation(self.q.conjugate())

    def display(self):
        print self.q.theta*2*180/math.pi
        print self.q.u.display()
