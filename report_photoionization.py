import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
import pandas as pd
from scipy.interpolate import CubicSpline
import warnings
import time
warnings.filterwarnings("ignore")
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.optimize import curve_fit
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import Pool
from tqdm import trange
num_cores = multiprocessing.cpu_count()
conv      = 27.2114 # One Hartree in electronvolts

def cross_section_photoionization(l,n,E): #in unit of pi*a0^2
        EK1 = E #entry in eV
        with open('Code/nh.txt', 'w') as f:
            f.write(str(n))
        with open('Code/lh.txt', 'w') as f:
            f.write(str(l))
        with open('Code/energy.txt', 'w') as f:
            f.write(str(EK1)+'d0')
        os.system('cd Code/ && gfortran work.f -std=legacy && ./a.out')
        with open('Code/output.txt', 'r') as f:
            for line in f:
                values = line.split()
        values = float(values[0])
        return values
   
color = ['','','blue','red','green','orange','purple','yellow']
l     =0
for n in [2,3,4]:
    photo = []
    Epp   = np.linspace(conv*0.51/n**2,conv*100/27.2114,51)
    for Ephoto in Epp:
        photo.append(cross_section_photoionization(l,n,Ephoto))
    photoionization = CubicSpline(Epp,photo)
    def fc_photoionization(E,n):
        if E < conv*0.51/n**2:
            return 0
        else:
            return photoionization(E)
    ydata = [fc_photoionization(E,n) for E in Epp]
    plt.plot(Epp,ydata,color=color[n],label='$n='+str(n)+'$,$l='+str(l)+'$',zorder=2)
'''
l=1
for n in [2,3,4]:
    photo = []
    Epp   = np.linspace(conv*0.51/n**2,conv*10/27.2114,51)
    for Ephoto in Epp:
        photo.append(cross_section_photoionization(l,n,Ephoto))
    photoionization = CubicSpline(Epp,photo)
    def fc_photoionization(E,n):
        if E < conv*0.51/n**2:
            return 0
        else:
            return photoionization(E)
    ydata = [fc_photoionization(E,n) for E in Epp]
    plt.plot(Epp,ydata,color=color[n],ls='--',label='$n='+str(n)+'$,$l='+str(l)+'$')
    
l=2
for n in [3,4]:
    photo = []
    Epp   = np.linspace(conv*0.51/n**2,conv*10/27.2114,51)
    for Ephoto in Epp:
        photo.append(cross_section_photoionization(l,n,Ephoto))
    photoionization = CubicSpline(Epp,photo)
    def fc_photoionization(E,n):
        if E < conv*0.51/n**2:
            return 0
        else:
            return photoionization(E)
    ydata = [fc_photoionization(E,n) for E in Epp]
    plt.plot(Epp,ydata,color=color[n],ls='-.',label='$n='+str(n)+'$,$l='+str(l)+'$')
    
l=3
for n in [4]:
    photo = []
    Epp   = np.linspace(conv*0.51/n**2,conv*10/27.2114,51)
    for Ephoto in Epp:
        photo.append(cross_section_photoionization(l,n,Ephoto))
    photoionization = CubicSpline(Epp,photo)
    def fc_photoionization(E,n):
        if E < conv*0.51/n**2:
            return 0
        else:
            return photoionization(E)
    ydata = [fc_photoionization(E,n) for E in Epp]
    plt.plot(Epp,ydata,color=color[n],ls=':',label='$n='+str(n)+'$,$l='+str(l)+'$')
'''

plt.legend()
plt.xscale('log')
plt.yscale('log')

plt.xlabel('$E$ [eV]')
plt.ylabel('$\sigma_{\mathrm{photoionization}}$ [$\pi a_0^2$]')
plt.grid(True, which="both", ls="-")
plt.xlim([conv*0.51/4**2,conv*10/27.2114])
plt.xlim(5e-1,1e2)
plt.ylim([1e-4,1e0])
plt.fill_between([1e-1,1e2],1e-4,1e0,color='white',alpha=1,zorder=0)
plt.tight_layout()
plt.savefig('report_photoionization.png',dpi=500,transparent=True)
