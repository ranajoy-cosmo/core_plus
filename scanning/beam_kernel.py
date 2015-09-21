import numpy as np
import healpy as hp

def gaussian_2d(x,y,ux,uy,fwhm_x,fwhm_y):
    factor = 2*np.sqrt(2*np.log(2))
    sigma_x = fwhm_x/factor
    sigma_y = fwhm_y/factor
    norm = 2*np.pi*sigma_x*sigma_y
    t1 = np.exp(-0.5*((x - ux)/sigma_x)**2)
    t2 = np.exp(-0.5*((y - uy)/sigma_y)**2)
    return (1/norm)*t1*t2

def get_mesh(nside, fwhm, cutoff = 4):
    factor = 2*np.sqrt(2*np.log(2))
    fwhm_max = np.max(fwhm)
    sigma = fwhm_max/factor
    dd = hp.nside2resol(nside, arcmin = True)
    nxy = int(cutoff*sigma/dd)
    x = np.arange(-1.0*nxy*dd, nxy*dd, dd)
    y = np.arange(-1.0*nxy*dd, nxy*dd, dd)
    return np.meshgrid(x,y)

del_beta = 1.71*(np.arange(9) - 4)
beam_kernel = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 1.0])
beam_kernel = beam_kernel/np.sum(beam_kernel)

