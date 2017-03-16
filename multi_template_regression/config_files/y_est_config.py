#import healpy as hp
#import os
#import numpy as np
#from simulation.lib.utilities.generic_class import Generic
#from simulation.global_config import global_paths
from default_config import *
#
#config = Generic()
#
#BOLO_CONFIG_FOLDER = "simulation.template_subtraction.bolo_config_files."
#
##*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
## Sim Action 
##*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.sim_tag = "bandpass_4_bolo_set_1"
#
##*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
## Read & Write
##*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#
config.bolo_list = ["bolo_0001"]
config.bolo_template_list = ["TEMPLATE_TD_0001", "TEMPLATE_SYNC_0001"]
#
#config.base_dir = global_paths.base_dir 
#config.general_data_dir = global_paths.output_dir
