#! /usr/bin/env python 

import numpy as np
import sys

def pad_array(a, pad_length, pad_value=0):
    b = np.full(a.size + 2*pad_length, pad_value)
    b[pad_length:-pad_length] = a
    return b

def convolve(a,b):
    len_a = a.size
    len_b = b.size
    if len_b%2 is 0:
        print "Length of b should be odd"
        sys.exit()
    pad_length = len_b/2
    a = pad_array(a, pad_length, 0)
    c = np.zeros(a.size)
    for i in range(-pad_length, pad_length + 1):
        c += np.roll(a, i)*b[pad_length + i]

    return c[pad_length:-pad_length]

