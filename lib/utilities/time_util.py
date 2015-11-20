#!/usr/bin/env python 

import time

def get_time_stamp():
    format = '%Y_%m_%d__%H_%M_%S'
    t = time.gmtime()
    return time.strftime(format, t)
