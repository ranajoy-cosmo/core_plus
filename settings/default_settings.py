import os
from simulation.lib.utilities.generic_class import Generic

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global path settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
global_paths = Generic()

#global_paths.base_dir = os.path.abspath("../")
global_paths.base_dir = "/home/ranajoy/COrE/simulation"

global_paths.output_dir = os.path.join(global_paths.base_dir, 'output')
global_paths.spectra_dir = os.path.join(global_paths.base_dir, 'spectra')
global_paths.maps_dir = os.path.join(global_paths.base_dir, 'maps', 'maps')

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global settings for scanning
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
global_scanning = Generic()

global_scanning.nside_in = 1024
global_scanning.nside_out = 256
global_scanning.lmax = 4000

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global system settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
global_system = Generic()

global_system.machine = 'Ranajoy_Macbook'
global_system.n_available_nodes = 1
global_system.n_cores_per_node = 4
global_system.n_cores_per_socket = 4
global_system.memory_per_node = 8                   #Giga-bytes
