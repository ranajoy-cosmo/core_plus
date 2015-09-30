import numpy as np
from simulation.lib.utilities.generic_class import Generic

settings = Generic()

settings.nside = 1024

settings.generate_pointing = False

base_folder = "/Users/banerji/CORE+/simulation/"

settings.data_folder = base_folder + "flight_data/segment_0001/"
settings.display_map = False 
settings.write_map = True
settings.pipe_with_simulation = True
