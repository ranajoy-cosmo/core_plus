#!/usr/bin/env python

from rotation_fast import *

f=open('viewing','w')
f_p=open('axis','w')
f_pol=open('pol_dir','w')

pi=math.pi
alpha=45*pi/180          #radians
beta=45*pi/180          #radians

theta_pol=90*pi/180

dT=1.0          #seconds
Tp=4*24*60*60.0          #seconds
Tr=60.0          #seconds

wp=2*math.pi/Tp          #radians/sec
wr=2*math.pi/Tr          #radians/sec

d_theta_p=wp*dT
d_theta_r=wr*dT

axis_prec=(1.0,0,0)          #precession about x-axis

axis_rot=(math.cos(alpha),0,math.sin(alpha))          #initial direction of rotation axis

v=(math.cos(alpha+beta),0,math.sin(alpha+beta))          #initial direction of pointing

for i in xrange(int(Tp/dT)+1):
    f_p.write(str(axis_rot[1])+'\t'+str(axis_rot[2])+'\n')
    #f.write(str(int(i/Tr)+1)+'\t'+str(i%int(Tr)+1)+'\t'+str(v[0])+'\t'+str(v[1])+'\t'+str(v[2])+'\n')
    f.write(str(v[0])+'\t'+str(v[1])+'\t'+str(v[2])+'\n')
    #pol_dir=rotate(normalise(cross(axis_rot,v)),theta_pol,v)
    #f_pol.write(str(int(i/Tr)+1)+'\t'+str(i%int(Tr)+1)+'\t'+str(pol_dir[0])+'\t'+str(pol_dir[1])+'\t'+str(pol_dir[2])+'\n')
    v=rotate(v,d_theta_r,axis_rot)
    axis_rot=rotate(axis_rot,d_theta_p,axis_prec)
    #axis_rot,v=rotate_rk2_m(axis_prec,axis_rot,v,wp,wr,dT)
    
f.close()
f_p.close()
f_pol.close()
