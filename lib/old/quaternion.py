#!/usr/bin/env python

import math
from vector import Vector

class Quaternion:
    def __init__(self,q0=None,q_v=None,norm=None,theta=None):
        if q0==None:
            self.norm=norm
            self.theta=theta*math.pi/180.0
            self.u=q_v.get_unit()
            self.switch_rep()
        else:
            self.theta=None
            self.q0=q0
            self.q=q_v
            self.norm_sq=q0*q0+q_v.x*q_v.x+q_v.y*q_v.y+q_v.z*q_v.z

    def switch_rep(self):
            self.q0=self.norm*math.cos(self.theta)
            self.q=self.u.scalar_mul(self.norm*math.sin(self.theta))
            self.norm_sq=self.norm*self.norm
            

    def __add__(a,b):
        return Quaternion(a.q0+b.q0,Vector(a.q.x+b.q.x,a.q.y+b.q.y,a.q.z+b.q.z))
    
    def __sub__(a,b):
        return Quaternion(a.q0-b.q0,Vector(a.q.x-b.q.x,a.q.y-b.q.y,a.q.z-b.q.z))

    def __mul__(a,b):
        res_q0=a.q0*b.q0-a.q.inner(b.q)
        res_q=b.q.scalar_mul(a.q0)+a.q.scalar_mul(b.q0)+a.q.cross(b.q)
        return Quaternion(q0=res_q0,q_v=res_q)

    def conjugate(self):
        return Quaternion(self.q0,self.q.scalar_mul(-1.0))

    def inverse(self):
        res_q0=self.q0/self.norm_sq
        res_q=self.q.scalar_mul(-1.0/self.norm_sq)
        return Quaternion(res_q0,res_q)

    def display(self):
        if self.theta is None:
            print '('+str(self.q0)+','+str(self.q.x)+','+str(self.q.y)+','+str(self.q.z)+')'
        else:
            print math.sqrt(self.norm_sq)
            self.q.get_unit().display()
            print math.acos(self.q0/math.sqrt(self.norm_sq))*180/math.pi


        
            
            

    
    
