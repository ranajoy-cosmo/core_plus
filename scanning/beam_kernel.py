import numpy as np

beam_kernel = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 1.0])
beam_kernel = beam_kernel/np.sum(beam_kernel)

