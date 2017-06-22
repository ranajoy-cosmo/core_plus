import numpy as np
from scipy.constants import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
#pylab.rcParams['xtick.major.pad']='10'
#pylab.rcParams['ytick.major.pad']='10'

def factor_dust(nu,T,nu0,factor):
    x = h*nu/(k*T)
    x0 = h*nu0/(k*T)
    factor_dust = (nu/nu0)**(3*factor)*(np.exp(x0) - 1)/(np.exp(x) - 1)
    return factor_dust

def convT2I(nu): # dB/dT
    T0 = 2.725
    x = h*nu/(k*T0)
    conv = 2*h**2*nu**4/k/T0**2*np.exp(x)/((c*(np.exp(x) - 1))**2)*1e20
    return conv

Ni = 370      
GHz = 1e9
beta = 1.7   # dust
T_dust = 18  # K
nu0 = 140*GHz
#=================================================

## ## define filters: g(nu)
N_nu = 15000    # full ranges
N_nu1 = 10000   # sample in Delta_nu relate to delt_nu
Delta_nu = 0.25*nu0
delt_nu = Delta_nu/N_nu1  # small delta

nu_i = np.arange(nu0-Delta_nu/2 - (N_nu-N_nu1)/2*delt_nu, nu0+Delta_nu/2 + (N_nu-N_nu1)/2*delt_nu,delt_nu)

g = np.zeros((N_nu,Ni))
I_dust_nu = np.zeros(Ni)
integral_cmb = np.zeros(Ni)
T_dust_factor = np.zeros(Ni)
factor_bandpass = np.zeros(Ni)
np.random.seed(Ni)
#fig = plt.figure(figsize=(6,4),dpi=100)

for i in range(Ni):
    #print 'detector band pass',i    
    e_nu1 = np.random.rand(1)*0.01*nu0
    e_nu2 = np.random.rand(1)*0.01*nu0   
    p1 = long((N_nu-N_nu1)/2 + (e_nu1)/delt_nu)
    p2 = long((N_nu-N_nu1)/2 + (Delta_nu + e_nu2)/delt_nu)    
    g[p1:p2,i] = 1
    
    for j in range(N_nu):
        I_dust_nu[i] += g[j,i]*factor_dust(nu_i[j],T_dust,nu0,beta)*delt_nu  / (nu_i[j]**2)
        integral_cmb[i] += g[j,i]*convT2I(nu_i[j])*delt_nu  / (nu_i[j]**2)       
    factor_bandpass[i] = I_dust_nu[i] / integral_cmb[i]
    T_dust_factor[i] = factor_bandpass[i]*convT2I(nu0)
print 'Dust factor of',Ni,'detectors is',T_dust_factor

#plt.plot(nu_i,g,linewidth=1.5)
#plt.legend(loc='lower center',fontsize=16)
#plt.title('boxcar bandpass mismatch', fontsize=27)
#plt.xlabel(r'$\nu$[Hz] ',fontsize=30)
#plt.tick_params(axis='x', labelsize=28, which='major', pad=10)
#plt.tick_params(axis='y', labelsize=28, which='major', pad=10)
#plt.tick_params(which='both', width=2)
#plt.tick_params(which='major', length=7)
#plt.tick_params(which='minor', length=7)
#plt.grid()
#plt.savefig('bandpas_function',dpi=400)
#plt.show()
print ' standard deviation of dust number',np.std(T_dust_factor),'and mean value', np.mean(T_dust_factor)





