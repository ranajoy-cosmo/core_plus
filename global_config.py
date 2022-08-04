import os
from simulation.lib.utilities.generic_class import Generic
# Testing line change

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global system settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_system = Generic()

hosts = ["edison", "cori", "apcmc", "apccl"] 
try:
    global_system.host = filter(lambda h : h in os.environ["HOST"], hosts)[0]
except KeyError:
    global_system.host = filter(lambda h : h in os.environ["HOSTNAME"], hosts)[0]

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Global path settings
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

global_paths = Generic()

home_dir = os.environ["HOME"]
base_dir_name = "simulation"

if global_system.host is "apcmc":
    global_paths.base_dir = os.path.join(home_dir, "CORE+/simulation")
    global_paths.scratch_dir = global_paths.base_dir
    global_paths.output_dir = os.path.join(global_paths.scratch_dir, 'output')
    global_paths.long_term_output_dir = os.path.join(global_paths.scratch_dir, 'long_term_output')
    global_paths.maps_dir = os.path.join(global_paths.scratch_dir, 'maps')

if global_system.host is "apccl":
    global_paths.base_dir = os.path.join(home_dir, base_dir_name) 
    global_paths.scratch_dir = "/workdir/banerji" 
    global_paths.output_dir = os.path.join(global_paths.scratch_dir, 'core_output')
    global_paths.long_term_output_dir = os.path.join(global_paths.scratch_dir, 'core_long_term_output')
    global_paths.maps_dir = os.path.join(global_paths.scratch_dir, 'core_maps')

if global_system.host in ["edison", "cori"]:
    global_paths.base_dir = os.path.join(home_dir, base_dir_name)
    global_paths.scratch_dir = os.environ["SCRATCH"]
    global_paths.general_output_dir = os.path.join(global_paths.scratch_dir, 'core_output')
    global_paths.long_term_output_dir = os.path.join(global_paths.scratch_dir, 'core_long_term_output')
    global_paths.maps_dir = os.path.join(global_paths.scratch_dir, 'core_maps')

global_paths.spectra_dir = os.path.join(global_paths.base_dir, 'spectra')
global_paths.simulation_code_dir = os.path.join(global_paths.base_dir, 'timestream_simulation')
global_paths.map_making_code_dir = os.path.join(global_paths.base_dir, 'map_maker')
