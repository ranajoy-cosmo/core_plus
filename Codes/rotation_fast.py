#!/usr/bin/env python

import math
from math import cos,sin

def add(a,b):
    return (a[0]+b[0],a[1]+b[1],a[2]+b[2])

def dot(a,b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def cross(a,b):
    return (a[1]*b[2]-b[1]*a[2],a[2]*b[0]-b[2]*a[0],a[0]*b[1]-b[0]*a[1])

def scalar_multiply(a,v):
    return (a*v[0],a*v[1],a*v[2])

def normalise(a):
    n=math.sqrt(dot(a,a))
    return (a[0]/n,a[1]/n,a[2]/n)

def rotate(v,theta,u):
    #if dot(u,u) is not 1:
     #   u=normalise(u)

    res=scalar_multiply(cos(theta),v)
    temp=(1-cos(theta))*dot(u,v)
    res=add(res,scalar_multiply(temp,u))
    return add(res,scalar_multiply(sin(theta),cross(u,v)))

def rotate_euler(v,theta,u):
    dv=scalar_multiply(theta,cross(u,v))
    res=add(v,dv)
    #return normalise(res)
    return res

def f(w,u,v):
    return scalar_multiply(w,cross(u,v))

def rotate_rk2(v,w,dt,u):
    f0=f_v(w,u,v)
    f1=f_v(w,u,add(v,scalar_multiply(dt/2,f0)))
    f2=f_v(w,u,add(v,scalar_multiply(dt/2,f1)))
    f3=f_v(w,u,add(v,scalar_multiply(dt,f2)))
    
    dv=add(add(f0,scalar_multiply(2,f1)),add(scalar_multiply(2,f2),f3))
    
    return add(v,scalar_multiply(dt/6,dv))

def rotate_rk2_m(p,r,v,w_p,w_r,dt):
    a_r=f(w_p,p,r)
    b_r=f(w_p,p,add(r,scalar_multiply(dt/2,a_r)))
    c_r=f(w_p,p,add(r,scalar_multiply(dt/2,b_r)))
    d_r=f(w_p,p,add(r,scalar_multiply(dt,c_r)))
    
    a_v=f(w_r,r,v)
    b_v=f(w_r,add(r,scalar_multiply(dt/2,a_r)),add(v,scalar_multiply(dt/2,a_v)))
    c_v=f(w_r,add(r,scalar_multiply(dt/2,b_r)),add(v,scalar_multiply(dt/2,b_v)))
    d_v=f(w_r,add(r,scalar_multiply(dt,c_r)),add(v,scalar_multiply(dt,c_v)))

    dr=add(add(a_r,scalar_multiply(2,b_r)),add(scalar_multiply(2,c_r),d_r))
    dv=add(add(a_v,scalar_multiply(2,b_v)),add(scalar_multiply(2,c_v),d_v))

    return add(r,scalar_multiply(dt/6,dr)),add(v,scalar_multiply(dt/6,dv))
    

    
