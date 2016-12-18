from default_config import *

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.simulate_ts = True/False. True means it will simulate the entire timestream signal and pointing from scratch
config.simulate_ts = False

config.scan_tag = "scan" 

#Options for input_pol_type = ["TQU", "QU", "T"]
#Default -> "TQU"
# A _ means that the input map has a component which will not be read"
config.pol_type = "TQU"
config.map_making_tag = "rec_pair_1"

config.take_diff_signal = False
config.subtract_template = False

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#If simulate_ts is true, this is used for simulation as well as map_making, otherwise these are the data segments that are read.

config.bolo_config_file = BOLO_CONFIG_FOLDER + "bolo_config" 
config.bolo_list = ['bolo_0001a', 'bolo_0001b']

config.segment_list = range(120)
