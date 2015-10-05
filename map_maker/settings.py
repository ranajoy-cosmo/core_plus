import numpy as np
from simulation.lib.utilities.generic_class import Generic

settings = Generic()

settings.nside = 1024

settings.generate_pointing = False

base_folder = "/Users/banerji/CORE+/simulation/output/"
tag = "test"
settings.output_folder = base_folder + tag + '/'
settings.display_map = False 
settings.write_map = True
settings.pipe_with_simulation = True
