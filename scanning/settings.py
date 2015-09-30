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
settings.plot_beam = True 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# I/O Params 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

base_folder = "/Users/banerji/CORE+/simulation/"
settings.input_map_path = base_folder + "maps_and_spectra/maps/test_map.fits"
settings.data_folder = base_folder + "flight_data/segment_0001/"

settings.write_signal = False 
settings.display_scanned_map = False
settings.write_scanned_map = True
settings.pipe_with_map_maker = True 
