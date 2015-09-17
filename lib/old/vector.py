#!/usr/bin/env python

import math

class Vector:
    def __init__(self,x=None,y=None,z=None):
        self.x=x
        self.y=y
        self.z=z
        if x is None: self.norm_sq=None
        else: self.norm_sq=x*x+y*y+z*z

    def get_data(self,x=0,y=0,z=0):
        self.x=x
        self.y=y
        self.z=z
        self.norm_sq=x*x+y*y+z*z

    def __add__(a,b):
        return Vector(a.x+b.x,a.y+b.y,a.z+b.z)

    def __sub__(a,b):
        return Vector(a.x-b.x,a.y-b.y,a.z-b.z)

    def norm(self):
        return math.sqrt(self.norm_sq)

    def scalar_mul(self,s):
        return Vector(s*self.x,s*self.y,s*self.z)
    
    def normalise(self):
        if self.norm_sq is not 1:
            norm=math.sqrt(self.norm_sq)
            self.x=self.x/norm
            self.y=self.y/norm
            self.z=self.z/norm
        else: pass

    def get_unit(self):
        if self.norm_sq is not 1:
            norm=math.sqrt(self.norm_sq)
            return Vector(self.x/norm,self.y/norm,self.z/norm)
        else: return self
    
    def inner(self,b):
        return self.x*b.x+self.y*b.y+self.z*b.z

    def cross(self,b):
        return Vector(self.y*b.z-b.y*self.z,self.z*b.x-b.z*self.x,self.x*b.y-b.x*self.y)

    def display(self):
        print '('+str(self.x)+','+str(self.y)+','+str(self.z)+')'


