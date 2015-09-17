#!/usr/bin/env python

from vector import Vector
from quaternion import Quaternion
from rotation import Rotation
import math

f=open('viewing','w')
f_p=open('axis','w')

pi=math.pi
alpha=45*pi/180          #radians
beta=45*pi/180          #radians

dT=1          #seconds
Tp=4*24*60*60*1.0         #seconds
Tr=60.0          #seconds

wp=360.0/Tp          #degrees/sec
wr=360.0/Tr          #degrees/sec

d_theta_p=wp*dT
d_theta_r=wr*dT

axis_prec=Vector(1.0,0,0)          #precession about x-axis

axis_rot=Vector(math.cos(alpha),0,math.sin(alpha))          #initial direction of rotation axis

v=Vector(math.cos(alpha+beta),0,math.sin(alpha+beta))          #initial direction of pointing


rot_p=Rotation(theta=d_theta_p/2,d=axis_prec,mode=1)

for i in range(int(Tp/dT+1)):
    f_p.write(str(axis_rot.x)+'\t'+str(axis_rot.y)+'\t'+str(axis_rot.z)+'\n')
    axis_rot=rot_p.rotate(axis_rot)
    rot_r=Rotation(theta=d_theta_r,d=axis_rot,mode=1)
    f.write(str(v.x)+'\t'+str(v.y)+'\t'+str(v.z)+'\n')
    v=rot_r.rotate(v)
    axis_rot=rot_p.rotate(axis_rot)

f.close()
f_p.close()

