
from simulation.timestream_simulation.config_files.default_config import config
from simulation.params.custom_params import global_paths, global_scanning, global_system

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Sim Action 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.sim_tag = "test_runs"  
config.scan_tag = "noise_grad_T_test" 

#Options for input_pol_type = ["TQU", "T_only", "noise_only"]
#Default -> "TQU"
config.sim_pol_type = "TQU"

config.add_noise = True

config.do_pencil_beam = True

config.gal_coords = False 

config.bolo_config_file = "100_bolos_optimal"  # rename by 100 bolos

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Scan 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.t_year = 366*24*60*60.0                    #seconds
config.t_prec = 4*24*60*60.0                     #seconds
config.t_spin = 120.0                             #seconds
config.sampling_rate = 113                          #Hz

config.alpha = 30.0                                #degrees
config.beta = 65.0                                #degrees

config.oversampling_rate = 1

config.nside_in = 1024

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Noise 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Options for noise_type = ["white", "1_over_f"]
#config.noise_type = "white"
#config.noise_sigma = 1200                           #\mu K 
#config.noise_f_knee = 0
#config.noise_alpha = 0

#Options for noise_type = ["white_noise", "1_over_f"]
config.noise_type = "white_noise"
 
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Data selection
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.bolo_list = ['bolo_0001a', 'bolo_0001b', 'bolo_0002a', 'bolo_0002b','bolo_0003a', 'bolo_0003b', 'bolo_0004a', 'bolo_0004b',
                    'bolo_0005a', 'bolo_0005b', 'bolo_0006a', 'bolo_0006b','bolo_0007a', 'bolo_0007b', 'bolo_0008a', 'bolo_0008b',
                    'bolo_0009a', 'bolo_0009b', 'bolo_0010a', 'bolo_0010b','bolo_0011a', 'bolo_0011b', 'bolo_0012a', 'bolo_0012b',
                    'bolo_0013a', 'bolo_0013b', 'bolo_0014a', 'bolo_0014b','bolo_0015a', 'bolo_0015b', 'bolo_0016a', 'bolo_0016b',
                    'bolo_0017a', 'bolo_0017b', 'bolo_0018a', 'bolo_00018b','bolo_0019a', 'bolo_0019b', 'bolo_0020a', 'bolo_0020b',
                    'bolo_0021a', 'bolo_0021b', 'bolo_0022a', 'bolo_0022b','bolo_0023a', 'bolo_0023b', 'bolo_0024a', 'bolo_0024b',
                    'bolo_0025a', 'bolo_0025b', 'bolo_0026a', 'bolo_0026b','bolo_0027a', 'bolo_0027b', 'bolo_0028a', 'bolo_0028b',
                    'bolo_0029a', 'bolo_0029b', 'bolo_0030a', 'bolo_0030b','bolo_0031a', 'bolo_0031b', 'bolo_0032a', 'bolo_0032b',
                    'bolo_0033a', 'bolo_0033b', 'bolo_0034a', 'bolo_0034b','bolo_0035a', 'bolo_0035b', 'bolo_0036a', 'bolo_0036b',
                    'bolo_0037a', 'bolo_0037b', 'bolo_0038a', 'bolo_00038b','bolo_0039a', 'bolo_0039b', 'bolo_0040a', 'bolo_0040b',
                    'bolo_0041a', 'bolo_0041b', 'bolo_0042a', 'bolo_0042b','bolo_0043a', 'bolo_0043b', 'bolo_0044a', 'bolo_0044b',
                    'bolo_0045a', 'bolo_0045b', 'bolo_0046a', 'bolo_0046b','bolo_0047a', 'bolo_0047b', 'bolo_0048a', 'bolo_0048b',
                    'bolo_0049a', 'bolo_0049b', 'bolo_0050a', 'bolo_0050b']


config.t_segment = 12*60*60.0                      #seconds
config.segment_list = range(4)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Read & Write
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

#Available options for timestream_data_products = ["timestream_data", "scanned_map"]
#"timestream_data" -> will write the simulated pointing, polarisation angle and the simulated signal in the scan data directory
#"scanned_map" -> will write the input hitmap and scanned map in the scan data directory
config.timestream_data_products = ["timestream_data", "hitmap", "noise"]

config.base_dir = global_paths.base_dir 
config.general_data_dir = global_paths.output_dir
config.write_beam = True
config.beam_file_name = "test_beam"

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
# Beam Parameters
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

config.simulate_beam = True
config.beam_cutoff = 4                                                             #fwhm
config.check_normalisation = True
config.display_beam_settings = True
