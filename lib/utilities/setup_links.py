#!/usr/bin/env python

import os
import sys
from simulation.params.custom_params import global_system, global_paths

scratch = os.environ["SCRATCH"]

def make_links():
    os.symlink(global_paths.output_dir, os.path.join(global_paths.base_dir, "output"))
    os.symlink(os.path.join(scratch, "core_long_term_output"), os.path.join(global_paths.base_dir, "long_term_output"))
    os.symlink(global_paths.maps_dir, os.path.join(global_paths.base_dir, "maps/maps"))

def remove_links():
    os.unlink(os.path.join(global_paths.base_dir, "output"))
    os.unlink(os.path.join(global_paths.base_dir, "long_term_output"))
    os.unlink(os.path.join(global_paths.base_dir, "maps/maps"))

if __name__=="__main__":
    try:
        action = sys.argv[1]
    except IndexError:
        action = None
    if action=="link":
        make_links()
    elif action=="unlink":
        remove_links()
    else:
        print "Input argument must be either \"link\" or \"unlink\""
