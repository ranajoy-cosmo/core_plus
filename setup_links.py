#!/usr/bin/env python

import os
import sys
from simulation.global_config import global_system, global_paths

def make_links():
    os.symlink(global_paths.output_dir, os.path.join(global_paths.base_dir, "output"))
    os.symlink(global_paths.long_term_output_dir, os.path.join(global_paths.base_dir, "long_term_output"))
    os.symlink(global_paths.maps_dir, os.path.join(global_paths.base_dir, "maps"))

def remove_links():
    os.unlink(os.path.join(global_paths.base_dir, "output"))
    os.unlink(os.path.join(global_paths.base_dir, "long_term_output"))
    os.unlink(os.path.join(global_paths.base_dir, "maps"))

if __name__=="__main__":
    try:
        action = sys.argv[1]
        if action=="link":
            make_links()
        elif action=="unlink":
            remove_links()
    except IndexError:
        print "Input argument must be either \"link\" or \"unlink\""
