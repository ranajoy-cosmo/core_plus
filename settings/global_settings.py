import os
from simulation.lib.utilities.generic_class import Generic

global_paths = Generic()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global path settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_paths.base_folder = os.path.abspath("../")

global_paths.output_folder = os.path.join(global_paths.base_folder, "output")
global_paths.spectra_folder = os.path.join(global_paths.base_folder, "maps_and_spectra/spectra")
global_paths.maps_folder = os.path.join(global_paths.base_folder, "maps_and_spectra/maps")

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global settings for resolution
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#global_resolution.nside = 2048
