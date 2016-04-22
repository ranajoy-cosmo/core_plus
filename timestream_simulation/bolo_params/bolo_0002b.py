import numpy as np
from simulation.lib.utilities import generic_class

bolo_params = generic_class.Generic()

bolo_params.bolo_name = "bolo_0002b"

#The offset of the bolo from the centre of the F-O-V axis
bolo_params.del_x = 0.0         #arc-min
bolo_params.del_y = 0.0         #arc-min

#The angle the polarisation direction of the bolometer with the W-E axis measured anti-clockwise
bolo_params.pol_ang = 135.0       #degrees

#The fwhm of the beam along two axes
bolo_params.fwhm = 8.0        #arcmin
#bolo_params.fwhm = np.sqrt(8.0**2 - 6.0**2)
bolo_params.conv_fwhm = 2.5946893621894396
#bolo_params.conv_fwhm = 0.0

#The angle the beam major axis of the bolometer with the W-E axis measured anti-clockwise
bolo_params.beam_angle = 0.0    #degrees