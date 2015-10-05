#! /usr/bin/env python 

import numpy as np

beam = np.arange(10, dtype = np.float_)

signal = np.arange(100)

signal_convolved = np.convolve(signal, beam, mode = 'same')

beam_weight = np.sum(beam)

signal_convolved_normalised = signal_convolved/beam_weight

