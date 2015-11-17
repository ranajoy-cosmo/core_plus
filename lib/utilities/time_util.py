#!/usr/bin/env python 

import time

def get_time_stamp():
    format = '%H_%M_%S__%d_%m_%Y'
    t = time.gmtime()
    return time.strftime(format, t)
