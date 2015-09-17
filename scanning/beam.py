#! /usr/bin/env python

import numpy as np
import pysimulators as ps

XSIZE = 256
YSIZE = 256

FWHMX = 3
FWHMY = 6

beam = ps.gaussian(shape = (YSIZE,XSIZE), fwhm = (FWHMX,FWHMY))

f_beam = np.fft.fft2(beam)

f_beam = np.imag(f_beam)

tmp_beam_1 = np.concatenate((f_beam[YSIZE/2:-1,XSIZE/2:-1],f_beam[YSIZE/2:-1,0:XSIZE/2]), axis = -1)

tmp_beam_2 = np.concatenate((f_beam[0:YSIZE/2,XSIZE/2:-1],f_beam[0:YSIZE/2,0:XSIZE/2]), axis = -1)

ref_beam = np.concatenate((tmp_beam_1,tmp_beam_2))
