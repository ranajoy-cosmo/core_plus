import numpy as np
import healpy as hp
from simulation.lib.utilities import generic_class
from simulation.settings.custom_settings import global_scanning, global_paths

settings = generic_class.Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Resolution settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.nside = global_scanning.nside_in
settings.beam_resolution = hp.nside2resol(settings.nside, arcmin=True)            #arcmins
settings.beam_cutoff = 3                    #FWHMs

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam shape settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.do_pencil_beam = False
settings.fwhm = 8.0                   #arcmins
settings.scan_radius = 45.0                 #degrees

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Display settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
settings.plot_beam = True
settings.display_beam_settings = True
