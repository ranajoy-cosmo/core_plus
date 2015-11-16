import os
from simulation.lib.utilities.generic_class import Generic

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global path settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
global_paths = Generic()

global_paths.tag = 'test'
global_paths.base_dir = os.path.abspath("../")

global_paths.output_dir = os.path.join(global_paths.base_dir, 'output')
global_paths.spectra_dir = os.path.join(global_paths.base_dir, 'spectra')
global_paths.maps_dir = os.path.join(global_paths.base_dir, 'maps', 'maps')

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global settings for scanning
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
global_scanning = Generic()

global_scanning.nside_in = 4096
global_scanning.nside_out = 1024
global_scanning.lmax = 4000
