import numpy as np
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.global_settings import global_settings

settings = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan Parameters 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

settings.nside = 2048 
settings.load_pointing = False
settings.generate_pointing = True

settings.do_beam_kernel = True 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read/Write Settings 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

tag = "test"
base_folder = global_settings.base_folder 
input_map_folder = base_folder + "maps_and_spectra/maps/"
settings.input_map = input_map_folder + "test_map_2048_0.fits"

settings.output_folder = global_settings.output_folder 

settings.write_signal = True 
settings.display_scanned_map = True 
settings.write_scanned_map = True
settings.pipe_with_map_maker = True 
