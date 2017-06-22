# Configure Matplotlib options
from setup_matplotlib import *
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import math as m
import numpy

# Load fiducials
fiducial1 = np.loadtxt("LowEll_best_fit_2016.dat") # l, TT, EE, BB, TE, TB
fiducial2 = np.loadtxt("LowEll_lensedtotCls+r01.dat") # l, TT, EE, BB, TE, TB
fiducial3 = np.loadtxt("LowEll_tensCls+r01.dat") # l, TT, EE, BB, TE, TB
fiducial4 = np.loadtxt("LowEll_lensedtotCls+r001.dat") # l, TT, EE, BB, TE, TB
fiducial5 = np.loadtxt("LowEll_tensCls+r001.dat") # l, TT, EE, BB, TE, TB
fiducial6 = np.loadtxt("LowEll_lensedtotCls+r0001.dat") # l, TT, EE, BB, TE, TB
fiducial7 = np.loadtxt("LowEll_tensCls+r0001.dat") # l, TT, EE, BB, TE, TB

#Load data
Spettri_simulatiEE = np.loadtxt("Spettri_simulati_EE.dat") # l, EE, errEE




# Create the plot
width=16
fig = plt.figure(figsize=(cm2inch(width), cm2inch(width*6/8.)))
plt.subplots_adjust(hspace=0.001)
plt.subplots_adjust(wspace=0.001)

ax1 = fig.add_subplot(111)

#Plotting error bars (I plot them at first to not hide all the other lines)
plt.fill_between(Spettri_simulatiEE[:,0],Spettri_simulatiEE[:,1]-Spettri_simulatiEE[:,2],Spettri_simulatiEE[:,1]+Spettri_simulatiEE[:,2],facecolor='b',alpha=0.3,linewidth=0.0)

#plotting EE fiducial
plt.plot(fiducial1[0:5000,0], fiducial1[0:5000,2], 'k')

#plotting BB fiducial
plt.plot(fiducial1[0:3900,0], fiducial1[0:3900,3], 'k')
plt.plot(fiducial2[0:3900,0], fiducial2[0:3900,3], 'k')
plt.plot(fiducial3[0:3900,0], fiducial3[0:3900,3], 'k--')
plt.plot(fiducial4[0:3900,0], fiducial4[0:3900,3], 'k')
plt.plot(fiducial5[0:3900,0], fiducial5[0:3900,3], 'k--')
plt.plot(fiducial6[0:3900,0], fiducial6[0:3900,3], 'k')
plt.plot(fiducial7[0:3900,0], fiducial7[0:3900,3], 'k--')

#Plotting data results
plt.plot(Spettri_simulatiEE[:,0], Spettri_simulatiEE[:,1], 'r',label='EE')

# axes limits
plt.ylim([2.0*10.0**(-7), 10.0**2]);
plt.xlim([1.9, 2000]);
ax1.set_yscale('log')
ax1.set_xscale('log')
plt.ylabel(r"$\ell (\ell +1) C_{\ell} /2\pi [\mu K^2]$");
plt.xlabel(r"Multipole $\ell$");
ax1.yaxis.labelpad = 10*width/17.;
ax1.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

plt.text(2,0.0025,'r=0.1', color='k', fontsize=12)
plt.text(2,0.00025,'r=0.01', color='k', fontsize=12)
plt.text(2,0.000025,'r=0.001', color='k', fontsize=12)

# set vertical y axis ticklables
for ticklabel in ax1.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")

for item in ([ax1.yaxis.label] + [ax1.xaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(12)

# reduce white space around figure
plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

# save to pdf with right bounding box
plt.savefig("Test.pdf", bbox_inches='tight')
