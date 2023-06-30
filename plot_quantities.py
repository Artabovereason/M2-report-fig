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
num_cores            = multiprocessing.cpu_count()
nlist = [2,3,4,5]
color      = ['','','blue','red','green','orange','purple','yellow']
true_conv = np.pi*2.80241e-1/0.874225

Ec = ['','',1.9191032244366055, 0.7492480995322649, 0.39206881586005515, 0.23780428489604535, 0.15993395358220985, 0.11269926439307765, 0.08374981501026807, 0.06395190877740402, 0.049905999130713925]
for i in nlist:
    plt.plot([Ec[i],Ec[i]],[0,100],color=color[i], linestyle = ':',alpha=1)
    
for n in nlist:
    data  = pd.read_csv('output/sigma_n='+str(n)+'.txt',delimiter=' ',header=None)
    xdata = data[0]
    ydata = data[1]
    #plt.plot(xdata,[true_conv*y/n for y in ydata],color=color[n],label='$n='+str(n)+'$')

for n in nlist:
    data0  = pd.read_csv('output/inf_sigma_n='+str(n)+'_l=0.txt',delimiter=' ',header=None)
    data1  = pd.read_csv('output/inf_sigma_n='+str(n)+'_l=1.txt',delimiter=' ',header=None)
    xdata0 = data0[0]
    ydata0 = data0[1]

    xdata1 = data1[0]
    ydata1 = data1[1]
    plt.plot(xdata0,[true_conv*y/n for y in ydata0],ls='-',color=color[n],zorder=10,label='$n='+str(n)+'$')
    #if n != 5:
    #    plt.plot(xdata1,[true_conv*y/n for y in ydata1],ls=':',color=color[n])
    #plt.plot(xdata1,[true_conv*y/n for y in ydata1],ls='-.',color=color[n])

for n in nlist:
    data  = pd.read_csv('output/alt_sigma_n='+str(n)+'_l=0.txt',delimiter=' ',header=None)
    xdata = data[0]
    ydata = data[1]
    #plt.plot(xdata,[true_conv*y/n for y in ydata],color=color[n],label='$n='+str(n)+'$',alpha=0.3)


array = pd.read_csv('H(1s)-H(2s).csv')
xdata = [float(i) for i in array['xdata']]
ydata = [float(i) for i in array['ydata']]
plt.plot(xdata,ydata,'x',color='blue',label='H(1s)+H(2s)',zorder=0,alpha=1)

array = pd.read_csv('H(1s)-H(3s).csv')
xdata = [float(i) for i in array['xdata']]
ydata = [float(i) for i in array['ydata']]
plt.plot(xdata,ydata,'x',color='red',label='H(1s)+H(3s)',alpha=1)

array = pd.read_csv('H(1s)-H(4s).csv')
xdata = [float(i) for i in array['xdata']]
ydata = [float(i) for i in array['ydata']]
plt.plot(xdata,ydata,'x',color='green',label='H(1s)+H(4s)',alpha=1)
'''
array = pd.read_csv('n_10.csv')
xdata = [float(i)**2 for i in array['xdata']]
ydata = [float(i) for i in array['ydata']]
plt.plot(xdata,ydata,'-',marker='.',color='grey')

array = pd.read_csv('n_15.csv')
xdata = [float(i)**2 for i in array['xdata']]
ydata = [float(i) for i in array['ydata']]
plt.plot(xdata,ydata,'-',marker='.',color='grey')

array = pd.read_csv('n_20.csv')
xdata = [float(i)**2 for i in array['xdata']]
ydata = [float(i) for i in array['ydata']]
plt.plot(xdata,ydata,'-',marker='.',color='grey')

array = pd.read_csv('n_25.csv')
xdata = [float(i)**2 for i in array['xdata']]
ydata = [float(i) for i in array['ydata']]
plt.plot(xdata,ydata,'-',marker='.',color='grey')
'''

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E$ [eV]')
plt.ylabel('$\\sigma(E)$ [$10^{-16}$ cm$^2$]')
plt.xlim(1e-1,1e1)
plt.ylim(1e-2,100)
plt.legend()
plt.tight_layout()
plt.savefig('output/sigma.png',dpi=300)