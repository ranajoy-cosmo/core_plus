#! /usr/bin/env python 

import numpy as np
from pyoperators import DegreesOperator, RadiansOperator, HomothetyOperator

#Radians to Degrees
rad2deg = DegreesOperator()

#Radians to Arcminutes
rad2am = rad2deg*60.0

#Degrees to Radians
deg2rad = rad2deg.I

#Arcminutes to Radians
am2rad = rad2am.I

#Degrees to Arcminutes
deg2am = HomothetyOperator(60.0)

#Arcminutes to Degrees
am2deg = deg2am.I
