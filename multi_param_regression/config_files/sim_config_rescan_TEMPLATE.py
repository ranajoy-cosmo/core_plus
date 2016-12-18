from default_config import *

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
config.simulate_ts = True
config.pipe_with_map_maker = False

config.scan_tag = "scan"

#Options for input_pol_type = ["TQU", "QU", "T", "_QU", "noise_only"]
#Default -> "TQU"
# A _ means that the input map has a component which will not be read"
config.sim_pol_type = "QU"

config.bolo_config_file = BOLO_CONFIG_FOLDER + "bolo_config" 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.bolo_list = ['bolo_TEMPLATE_QU_re']
config.t_segment = 72*60*60.0                      #seconds
config.segment_list = range(120)
