#!/usr/bin/env python

import numpy as np
import healpy as hp
import sys
import simulation.power_spectrum.noise_power_spectra as nps
from simulation.lib.utilities.prompter import prompt
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

hitmap_file = "/global/homes/b/banerji/simulation/output/focal_plane_noise/36_bolo_4_year/hitmap.fits"
det_sens = 39.9
scan_freq = 85.0
noise_rms = det_sens * np.sqrt(scan_freq)
nside = 1024
lmax = 2000
fwhm = np.radians(7.7/60.0)
num_iteration_total = int(sys.argv[1])
num_iteration_per_rank = num_iteration_total / size 

hitmap = hp.read_map(hitmap_file)

sum_local = np.zeros(lmax+1)
sum_squared_local = np.zeros(lmax+1)

for i in xrange(num_iteration_per_rank):
    prompt("Rank : {}, doing iteration {} of {}\n".format(rank, i+1, num_iteration_per_rank))
    noise_map = np.random.normal(0.0, noise_rms, 12*1024**2) / np.sqrt(hitmap) / 2.0
    noise_spectra = hp.anafast(noise_map, lmax=lmax)
    sum_local += noise_spectra
    sum_squared_local += noise_spectra**2

sum = np.zeros(lmax+1)
sum_squared = np.zeros(lmax+1)

comm.Reduce(sum_local, sum, MPI.SUM, 0)
comm.Reduce(sum_squared_local, sum_squared, MPI.SUM, 0)

if rank == 0:
    spectra_mean = sum / num_iteration_total
    spectra_std = np.sqrt(sum_squared / num_iteration_total - spectra_mean**2)
    Bl = hp.gauss_beam(fwhm=fwhm, lmax=lmax, pol=False)
    spectra_mean /= Bl**2 
    spectra_std /= Bl**2
    spectra_th = nps.sensitivity_to_noise_mean_Cl(det_sens, fwhm, lmax)[0]
    np.save("noise_mc/spectra_mean", spectra_mean)
    np.save("noise_mc/spectra_sd", spectra_std)
    np.save("noise_mc/spectra_th", spectra_th)
