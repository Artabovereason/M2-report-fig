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

fine_structure_inv = 137.035999084
fine_structure     = 1/fine_structure_inv

### Allow to rasterize the data of the .xlsx in the correct python format
def data_to_data(array_data):
    ccache=[]
    for i in range(len(array_data)):
        if array_data[i] == array_data[i]:
            cache = ''
            for j in range(len(array_data[i])):
                if array_data[i][j] == '+' or array_data[i][j] == '-' and j!=0:
                    cache += 'E'
                cache += array_data[i][j]
            ccache.append(float("{:.8f}".format(float(cache))))
    return ccache


def R_trans(R,array_data):
    R_cache = []
    for i in range(len(array_data)):
        R_cache.append(R[i])
    return R_cache

n1               = 1 # lower state                               # Principal quantum number of the transition
df1              = pd.read_excel(os.path.join("dipole_strengths_data.xlsx"),engine='openpyxl',)
D                = data_to_data(np.array(df1['(1G,1U)'])) # Dipole moment Q in units of [e*a0] where e is the electronic charge
RR               = R_trans(np.array(df1['R']),D)          # Internuclear separation in bohrs, [a0]
dipolar_coupling = CubicSpline(RR, D)                   # Spline of the dipolar coupling of the (u,g) transition symmetry

'''
dipolar_coupling_tot : return the  absolute value of the total dipolar coupling between the two states we are considering.
Unit is [e*a0].
'''

def dipolar_coupling_tot(R):
    return dipolar_coupling(R)

conv     = 27.2114 # One Hartree in electronvolts
'''
U1        : return the potential of the lower state [eV], to use Eh, put conv = 1.
U2        : return the potential of the uper state [eV], to use Eh, put conv = 1.
U12       : return the difference of potential between the two states [eV], to use Eh, put conv = 1.
W         : ionization rate calculation
vrad      : radial collision velocity
first_int : return either Eq. 2.11 or Eq. 2.12 : the ionization probability
'''

if __name__ == '__main__':
    start_time = time.time()
    energies   = np.logspace(np.log(1e-2),np.log(10),301,base=np.exp(1))
    mu         = 1836.152/2.0 
    color      = ['','','blue','red','green','orange','purple','yellow']
    colorb     = ['','','cyan','pink','olive']
    R_plot     = np.linspace(1,16,1001)
    D12 = []
    for R in R_plot:
        D12.append(dipolar_coupling_tot(R)**2)
    plt.plot(R_plot,D12,color='black')
    #plt.legend(loc='upper left')
    #plt.xscale('log')
    plt.xlabel('$R$ [$a_0$]')
    plt.ylabel('$D_{12}^2(R)$ [$e^2a_0^2$]')
    plt.xlim(min(energies),max(energies))
    plt.ylim(0,1.1*max(D12))
    plt.xlim(0,16)
    plt.fill_between([0,17],0,1.1*max(D12),color='white',alpha=1)
    plt.tight_layout()
    plt.savefig('report_dipolar_coupling.png',dpi=500,transparent=True)
    print("--- %s seconds ---" % (time.time() - start_time))