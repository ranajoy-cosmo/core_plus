from simulation.lib.utilities import generic_class

bolo_params = generic_class.Generic()

bolo_params.bolo_name = "bolo_0003"

#The offset of the bolo from the centre of the F-O-V axis
bolo_params.del_x = 0.0         #arc-min
bolo_params.del_y = 0.0         #arc-min

#The angle the polarisation direction of the bolometer with the W-E axis measured anti-clockwise
bolo_params.pol_ang = 90.0       #degrees

#The fwhm of the beam along two axes
bolo_params.fwhm = 8.0        #arcmin
#bolo_params.ellipticity = 0.0
bolo_params.ellipticity = 0.05

#The angle the beam major axis of the bolometer with the W-E axis measured anti-clockwise
bolo_params.beam_angle = 0.0    #degrees
