from default_config import *

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#config.simulate_ts = True/False. True means it will simulate the entire timestream signal and pointing from scratch
config.simulate_ts = False

config.map_making_tag = "tm_grad_co_co"

#pol_type can be "TQU" OR "QU" OR "T"
config.pol_type = "QU"   

config.take_diff_signal = False
config.subtract_template = False
config.noise_only_map = False
 
config.notes = "Simulating the gradient along the scan direction and making the QU map from it"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#If simulate_ts is true, this is used for simulation as well as map_making, otherwise these are the data segments that are read.

config.bolo_list = ['bolo_0001a']

config.sim_type = "template"
config.template_type = "tm_gradient"
config.tm_gradient_type = ["tm_grad_co_co"]
