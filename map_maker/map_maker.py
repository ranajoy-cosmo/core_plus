#! /usr/bin/env python 

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from settings import settings
from pysimulators import ProjectionOperator
from pysimulators.sparse import FSRMatrix
import simulation.pointing.generate_pointing as gen_p

def make_map_from_signal():
    if settings.pipe_with_simulation:
        import simulation.scanning.scanning_T as scanning
        signal, P = scanning.simulate_beam_tod()

    sky_map = (P.T*P).I*P.T*signal

    if settings.display_map:
        hp.mollzoom(sky_map)
        plt.show()

    if settings.write_map:
        hp.write_map(settings.output_folder + "reconstructed_map.fits", sky_map)


if __name__=="__main__":
    make_map_from_signal()
    
