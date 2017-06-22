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
#
##*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
## Read & Write
##*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#
config.bolo_list = ["bolo_0001"]

#config.TEMPLATE_list = ["tm_grad_co", "tm_grad_cross", "tm_grad_co_co", "tm_grad_cross_cross", "tm_grad_co_cross"]
config.TEMPLATE_list = ["tm_grad_co_co", "tm_grad_cross_cross", "tm_grad_co_cross"]
#
#config.base_dir = global_paths.base_dir 
#config.general_data_dir = global_paths.output_dir
