import numpy as np
from scipy import fftpack
import pysimulators as ps
import matplotlib.pyplot as plt
from simulation.lib.plotting.my_imshow import new_imshow
import scipy as sp

#cmb_img = np.load("cmb_img.npy")

beam_1 = ps.datautils.gaussian(shape = (37, 37), fwhm = (8.0, 8.0))
beam_1 = beam_1/np.max(beam_1)
beam_2 = ps.datautils.gaussian(shape = (37, 37), fwhm = (8.0, 8.01))
beam_2 = beam_2/np.max(beam_2)

fft_beam_1 = fftpack.fft2(beam_1)
fft_beam_1[fft_beam_1 == 0] = 1e-18
fft_beam_2 = fftpack.fft2(beam_2)
fft_beam_2[fft_beam_2 == 0] = 1e-18

def draw_subplot(fig, plot_shape, plot_num, plot_object, title = None):
    ax = fig.add_subplot(plot_shape[0], plot_shape[1], plot_num)
    im = new_imshow(ax,plot_object)
    fig.colorbar(im, ax=ax)

fig = plt.figure()

draw_subplot(fig, (2,2), 1, beam_1, title = "Beam 1 : Symmetric")
draw_subplot(fig, (2,2), 2, beam_2, title = "Beam 2 : Asymmetric")
draw_subplot(fig, (2,2), 3, fftpack.fftshift(np.absolute(fft_beam_1)), title = "Beam 1 : FFT Absolute")
draw_subplot(fig, (2,2), 4, fftpack.fftshift(np.absolute(fft_beam_2)), title = "Beam 2 : FFT Absolute")

fig = plt.figure()

with np.errstate(divide='ignore', invalid='ignore'):
    f_conv_1 = fftpack.ifft2(fft_beam_1/fft_beam_1)
    f_conv_1_view = fftpack.fftshift(np.absolute(f_conv_1))
draw_subplot(fig, (1,2), 1, f_conv_1_view, title = "Convolution function 1/1")

with np.errstate(divide='ignore', invalid='ignore'):
    f_conv = fftpack.ifft2(fft_beam_1/fft_beam_2)
    f_conv_view = fftpack.fftshift(np.absolute(f_conv))
draw_subplot(fig, (1,2), 2, f_conv_view, title = "Convolution function 1/2")

"""
fig = plt.figure()

draw_subplot(fig, (1,1), 1, cmb_img)

img_conv_1 = sp.signal.convolve2d(cmb_img, beam_1, mode='same')/np.sum(beam_1)

fig = plt.figure()

draw_subplot(fig, (1,1), 1, img_conv_1)

img_conv_2 = sp.signal.convolve2d(cmb_img, beam_2, mode='same')/np.sum(beam_2)

fig = plt.figure()

draw_subplot(fig, (1,1), 1, img_conv_2)

img_conv_3 = sp.signal.convolve2d(cmb_img, np.array([[1]]), mode='same')

fig = plt.figure()

draw_subplot(fig, (1,1), 1, img_conv_3)


plt.show()
"""
