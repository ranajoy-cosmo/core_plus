import numpy as np
from simulation.lib.utilities.generic_class import Generic
from simulation.settings.global_settings import global_paths, global_scanning
import simulation.bolo.local_settings as scan_settings

settings = Generic()

settings.nside = global_scanning.nside

base_folder = global_paths.base_folder
settings.output_folder = global_paths.output_folder
settings.display_map = False 
settings.write_map = True

scan_params = scan_settings.settings
scan_params.write_signal = False
scan_params.pipe_with_map_maker = True
scan_params.display_scanned_map = False
