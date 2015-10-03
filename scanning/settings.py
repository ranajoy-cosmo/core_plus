import numpy as np
from simulation.lib.utilities.generic_class import Generic

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Parameters 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

settings.nside = 2048 
settings.load_pointing = False
settings.generate_pointing = True

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

settings.beam_cutoff = 4 
settings.normalise_beam = True 
settings.beam_fwhm = (8.0, 8.0)

settings.do_beam_kernel = True 
settings.plot_beam = False 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

tag = "test"
base_folder = "/Users/banerji/CORE+/simulation/"
input_map_folder = base_folder + "maps_and_spectra/maps/"
settings.input_map = input_map_folder + "test_map.fits"

settings.output_folder = base_folder + "output/" + tag + "/"

settings.write_signal = False 
settings.display_scanned_map = True 
settings.write_scanned_map = True
settings.pipe_with_map_maker = True 
