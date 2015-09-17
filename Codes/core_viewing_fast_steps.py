#!/usr/bin/env python

from rotation_fast import *

f=open('viewing','w')
f_p=open('axis','w')
f_pol=open('pol_dir','w')

pi=math.pi
alpha=45*pi/180          #radians
beta=45*pi/180          #radians

theta_pix=90*pi/180

dT=1.0          #seconds
Tp=4*24*60*60.0          #seconds
Tr=60.0          #seconds

wp=2*pi/Tp          #degrees/sec
wr=2*pi/Tr          #degrees/sec

d_theta_p=wp*Tr
d_theta_r=wr*dT

axis_prec=(1.0,0,0)          #precession about x-axis

axis_rot=(math.cos(alpha),0,math.sin(alpha))          #initial direction of rotation axis

v=(math.cos(alpha+beta),0,math.sin(alpha+beta))          #initial direction of pointing

for i in xrange(int(Tp/Tr)+1):
    v_prev=v
    f_p.write(str(axis_rot[0])+'\t'+str(axis_rot[1])+'\t'+str(axis_rot[2])+'\n')
    
    for j in xrange(int(Tr/dT)+1):
        pol_axis=rotate(normalise(cross(axis_rot,v)),theta_pix,v)
        #pol_axis=normalise(cross(axis_rot,v))
        f.write(str(v[0])+'\t'+str(v[1])+'\t'+str(v[2])+'\n')
        f_pol.write(str(i)+'\t'+str(j)+'\t'+str(pol_axis[0])+'\t'+str(pol_axis[1])+'\t'+str(pol_axis[2])+'\n')
        #v=rotate_rk2(v,wr,dT,axis_rot)
        v=rotate(v,d_theta_r,axis_rot)
    #axis_rot=rotate_rk2(axis_rot,wp,dT*Tr,axis_prec)
    axis_rot=rotate(axis_rot,d_theta_p,axis_prec)
    #v=rotate_rk2(v_prev,wp,dT*Tr,axis_prec)
    v=rotate(v_prev,d_theta_p,axis_prec)

f.close()
f_p.close()
f_pol.close()
