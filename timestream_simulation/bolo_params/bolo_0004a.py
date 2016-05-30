import numpy as np
from simulation.lib.utilities import generic_class

bolo_params = generic_class.Generic()

bolo_params.bolo_name = "bolo_0004a"

#The offset of the bolo from the centre of the F-O-V axis
bolo_params.del_x = 0.0         #arc-min
bolo_params.del_y = 0.0         #arc-min

#The angle the polarisation direction of the bolometer with the W-E axis measured anti-clockwise
bolo_params.pol_ang = 67.5       #degrees

#The fwhm of the beam along two axes
bolo_params.fwhm = 8.0        #arcmin
#bolo_params.conv_fwhm = 9.730085108210398                   #5% for 30'
bolo_params.conv_fwhm = 2.5946893621894396                 #5% for 8'
#bolo_params.conv_fwhm = 1.2973446810947198                  #5% for 4'
#bolo_params.conv_fwhm = 0.0

#The angle the beam major axis of the bolometer with the W-E axis measured anti-clockwise
bolo_params.beam_angle = 42.0    #degrees
