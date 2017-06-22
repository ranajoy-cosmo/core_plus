from default_config import *

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.simulate_ts = True/False. True means it will simulate the entire timestream signal and pointing from scratch
config.simulate_ts = False

config.map_making_tag = "rec_corrected"

config.sim_type = "signal"

#pol_type can be "TQU" OR "QU" OR "T"
config.pol_type = "QU"   

config.take_diff_signal = True 
config.subtract_template = True
config.noise_only_map = False

#config.TEMPLATE_list = ["tm_grad_co", "tm_grad_cross", "tm_grad_co_co", "tm_grad_cross_cross", "tm_grad_co_cross"]
config.TEMPLATE_list = ["tm_grad_co_co", "tm_grad_cross_cross", "tm_grad_co_cross"]
 
config.notes = "Making Q,U maps for the differenced signal of 0001a and 0001b"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#If simulate_ts is true, this is used for simulation as well as map_making, otherwise these are the data segments that are read.

config.bolo_list = ['bolo_0001']
