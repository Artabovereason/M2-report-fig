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
color      = ['','','blue','red','green','orange','purple','yellow','black','brown','pink']
from matplotlib.pyplot import figure
figure(figsize=(6, 6), dpi=500)
def photo_langevin(E,n):
    E = E/27.211
    return 2.80241e-1* (3*np.pi/(2**(2/3)))*((27*np.pi/16)*(1/E)*(1/(2+2*n**2))*(n**4)*(2*n**2+7))**(1/3)

energies = np.logspace(np.log(1e-2),np.log(10),31,base=np.exp(1))

for n in [2,3,4,5,6,7,8,9,10]:
    photo = [photo_langevin(E,n) for E in energies]
    plt.plot(energies,photo,color=color[n],label='$n='+str(n)+'$')

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

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E$ [eV]')
plt.ylabel('$\\sigma(E)$ [$10^{-16}$ cm$^2$]')
plt.xlim(1e-2,1e1)

plt.legend()
plt.tight_layout()
plt.savefig('output/sigma_langevin.png',dpi=300,transparent=True)
plt.clf()









