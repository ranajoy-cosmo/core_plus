import numpy as np
from simulation.lib.utilities.generic_class import Generic

settings = Generic()

settings.nside = 1024

settings.input_map_path = "/Users/banerji/CORE+/simulation/maps_and_spectra/maps/test_map.fits"

settings.generate_pointing = True

base_folder = "/Users/banerji/CORE+/simulation/"

settings.write_signal = False 
settings.data_folder = base_folder + "flight_data/segment_0001/"
settings.display_scanned_map = True
settings.write_scanned_map = True
settings.pipe_with_map_maker = True
