import numpy as np
from simulation.lib.utilities import generic_class
from simulation.timestream_simulation.custom_params import scan_params

bolo_params = generic_class.Generic()

bolo_params.bolo_name = "bolo_0001b"

#The offset of the bolo from the centre of the F-O-V axis
if scan_params.fixed_pointing_error == 0.0:
    bolo_params.pointing_offset_x = 0.0
    bolo_params.pointing_offset_y = 0.0
else:
    bolo_params.pointing_offset_x = np.random.normal(loc=0.0, scale=scan_params.fixed_pointing_error)          #arc-seconds
    bolo_params.pointing_offset_y = np.random.normal(loc=0.0, scale=scan_params.fixed_pointing_error)         #arc-seconds

#The angle the polarisation direction of the bolometer with the W-E axis measured anti-clockwise
bolo_params.pol_phase_ini = 90.0       #degrees

#The fwhm of the beam along two axes
bolo_params.fwhm = 8.0        #arcmin
#bolo_params.conv_fwhm = 9.730085108210398                   #5% for 30'
#bolo_params.conv_fwhm = 2.5946893621894396                 #5% for 8'
#bolo_params.conv_fwhm = 1.2973446810947198                  #5% for 4'
bolo_params.conv_fwhm = 0.0

#The angle the beam major axis of the bolometer with the W-E axis measured anti-clockwise
bolo_params.beam_angle = 5.0    #degrees
