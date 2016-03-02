#!/usr/bin/env python

import os
import sys
from simulation.params.custom_params import global_system, global_paths

scratch = os.environ["SCRATCH"]
os.symlink(global_paths.output_dir, os.path.join(global_paths.base_dir, "output"))
os.symlink(os.path.join(scratch, "core_long_term_output"), os.path.join(global_paths.base_dir, "long_term_output"))
os.symlink(global_paths.maps_dir, os.path.join(global_paths.base_dir, "maps/maps"))
