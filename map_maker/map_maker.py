#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix
import simulation.bolo.timestream_interp as ts
from local_settings import settings, scan_params
import os

def make_map_from_signal():

    signal, P = ts.do_simulation(scan_params)
    sky_map = (P.T*P).I*P.T*signal

    if settings.display_map:
        hp.mollzoom(sky_map)
        plt.show()

    if settings.write_map:
        hp.write_map(os.path.join(settings.output_folder, "reconstructed_map.fits"), sky_map)


if __name__=="__main__":
    make_map_from_signal()
    
