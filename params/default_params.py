import os
from simulation.lib.utilities.generic_class import Generic

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global system settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_system = Generic()

hosts = ["edison", "cori", "apc"] 
global_system.host = filter(lambda h : h in os.environ["HOST"], hosts)[0]

if global_system.host not in hosts:
    print "Host :", global_system.host, "not in list of known hosts. Please add in custom_settings"

if global_system.host is "apc":
    global_system.n_available_nodes = 1
    global_system.n_cores_per_node = 2
    global_system.memory_per_node = 8                   #Giga-bytes

if global_system.host is "edison":
    global_system.n_available_nodes = 1
    global_system.n_cores_per_node = 24
    global_system.n_cores_per_socket = 12
    global_system.memory_per_node = 64                   #Giga-bytes

if global_system.host is "cori":
    global_system.n_available_nodes = 1
    global_system.n_cores_per_node = 32
    global_system.n_cores_per_socket = 16
    global_system.memory_per_node = 128                   #Giga-bytes


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global path settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_paths = Generic()

home_dir = os.environ["HOME"]
base_dir_name = "simulation"

if global_system.host is "apc":
    global_paths.base_dir = os.path.join(home_dir, "CORE+/simulation")
    global_paths.output_dir = os.path.join(global_paths.base_dir, 'output')
    global_paths.maps_dir = os.path.join(global_paths.base_dir, 'maps')

if global_system.host in ["edison", "cori"]:
    global_paths.base_dir = os.path.join(home_dir, base_dir_name)
    scratch_dir = os.environ["SCRATCH"]
    global_paths.output_dir = os.path.join(scratch_dir, 'core_output')
    global_paths.maps_dir = os.path.join(scratch_dir, 'core_maps')

global_paths.spectra_dir = os.path.join(global_paths.base_dir, 'spectra')

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global settings for scanning
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
global_scanning = Generic()

global_scanning.nside_in = 4096
global_scanning.nside_out = 4096
global_scanning.lmax = 5000
