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
D               = data_to_data(np.array(df1['(1G,1U)'])) # Dipole moment Q in units of [e*a0] where e is the electronic charge
RR              = R_trans(np.array(df1['R']),D)          # Internuclear separation in bohrs, [a0]
dipolar_coupling = CubicSpline(RR, D)                   # Spline of the dipolar coupling of the (u,g) transition symmetry

'''
dipolar_coupling_tot : return the  absolute value of the total dipolar coupling between the two states we are considering.
Unit is [e*a0].
'''

def dipolar_coupling_tot(R):
    return dipolar_coupling(R)

conv     = 27.2114 # One Hartree in electronvolts

'''
df1      = pd.read_excel(os.path.join("adiabtic_curves.xlsx"),sheet_name='5p sigma u',engine='openpyxl')
cache    = df1['DISTANCE BOHR']
distance1 = [float(i) for i in cache]
cache    = df1['ENERGY HARTREE']
energy   = [conv*float(i) for i in cache]
U2_curve = CubicSpline(distance1, energy)                       # Adiabatic energy curve for the upper state [eV], to use Eh, put conv = 1.

df1      = pd.read_excel(os.path.join("adiabtic_curves.xlsx"),sheet_name='1s sigma g',engine='openpyxl')
cache    = df1['DISTANCE BOHR']
distance2= [float(i) for i in cache]
cache    = df1['ENERGY HARTREE']
energy   = [conv*float(i) for i in cache]
U1_curve = CubicSpline(distance2, energy)                     # Adiabatic energy curve for the lower state [eV], to use Eh, put conv = 1.
'''

array     = pd.read_csv('discurves/0_0_1_new.dat',delimiter='    ',engine='python')
xdata     = [float(i) for i in array['xdata']]
ydata     = [float(i)+0.5 for i in array['ydata']]
inv_xdata = [1/i for i in xdata]
ydata     = [conv*(ydata[i]+inv_xdata[i]) for i in range(len(ydata)) if xdata[i] > 0]
U1_curve  = CubicSpline(xdata, ydata)

array     = pd.read_csv('discurves/1_0_1.dat',delimiter='    ',engine='python')
xdata     = [float(i) for i in array['xdata']]
ydata     = [float(i)+0.5 for i in array['ydata']]
inv_xdata = [1/i for i in xdata]
ydata     = [conv*(ydata[i]+inv_xdata[i]) for i in range(len(ydata)) if xdata[i] > 0]
U2_curve  = CubicSpline(xdata, ydata)


'''
U1        : return the potential of the lower state [eV], to use Eh, put conv = 1.
U2        : return the potential of the uper state [eV], to use Eh, put conv = 1.
U12       : return the difference of potential between the two states [eV], to use Eh, put conv = 1.
W         : ionization rate calculation
vrad      : radial collision velocity
first_int : return either Eq. 2.11 or Eq. 2.12 : the ionization probability
'''

def U1(R):
    return  U1_curve(R) 

def U2(R):
    return  U2_curve(R) 

def U12(R):
    return U2(R)-U1(R)

# Find maximal R region
def find_RN(n):
    RN = 16
    while U12(RN)< conv*0.5/n**2:
        RN -= 0.01
    return RN

def vrad(rho,E,R,mu,shift):
        return np.sqrt(2/mu)*np.sqrt(np.abs(E*(1-(rho/R)**2)-(U2(R)-shift)))

# Find classical turning point
def find_RO(rho,E,shift):
    R0 = 1
    while (E-(U2(R0)-shift)-E*(rho/R0)**2) <0:
        R0 += 0.001
        if R0 > 16:
            break
    return R0

def find_RE(E):
    R = 16
    while U2(R)-U1(R) < E:
        R -= 0.01
        if R < 1:
            break
    return R

if __name__ == '__main__':
    start_time = time.time()
    energies   = np.logspace(np.log(5e-1),np.log(10),31,base=np.exp(1))
    Rcl        = ['','',3.709 ,4.789,5.527, 6.085, 6.526, 6.895, 7.21, 7.489, 7.732]
    color      = ['','','blue','red','green','orange','purple','yellow']
    colorb     = ['','','cyan','pink','olive']

    def cross_section_photoionization(l,n,E): #in unit of pi*a0^2
        EK1 = E#*27.2114#EK(R,n)*27.2114 #because the input format of the code is like that, should change that in the future to avoid confusion
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
   
    conv_ato_cm2 = 0.874225 
    mu           = 1836.152/2.0 #mp/2 in unit of me
    c            = 137.037 
    
    def job(E,n,RN):
        rho = 0
        #RN = find_RN(n)
        while find_RO(rho,E,0) < RN:
            rho += 0.01
            if rho > 16:
                break
        if rho == 0:
            return rho
        return rho

    color = ['','','blue','red','green','orange','purple','yellow']
    Ec    = ['','',2,6e-1,4e-1,2e-1]
    def fc_pci(E,n,RN,rho,RE):
        ccache = []
        cccache = []
        for R in np.linspace(find_RO(rho,E,0),RE,5001):
            top = c*(dipolar_coupling_tot(R)**2) * fc_photoionization(U12(R),n)*U12(R)**3/(2*np.pi)  
            #top = c*(dipolar_coupling_tot(R)**2) * fc_photoionization(U12(R),n)*E**3/(2*np.pi)  
            ccache.append(-2*(conv**(-3)/np.sqrt(conv**(-1)))*(top/vrad(rho,E,R,mu,0)))
        first  = 0.5*(1-np.exp(np.sum(ccache)/len(ccache)))#*rho

        for R in np.linspace(RE,RN,5001):
            top = c*(dipolar_coupling_tot(R)**2) * fc_photoionization(U12(R),n)*U12(R)**3/(2*np.pi)  
            cccache.append(-(conv**(-3)/np.sqrt(conv**(-1)))*(top/vrad(rho,E,R,mu,0)))
        second = np.exp(-np.sum(cccache)/len(cccache))
    
        return first*second*rho
    def fc_sigma(E,n,RN,RE):
        rhomaxE= job(E,n,RE)
        rhos   = np.linspace(0,rhomaxE,101)
        '''
        rhomaxE= job(E,n,RE)
        if RN < RE:
            rhos   = np.linspace(0,rhomax,101)
        else:
            rhos   = np.linspace(0,rhomaxE,101)
        '''
        cache  = []
        for rho in rhos:
            cache.append(fc_pci(E,n,RN,rho,RE))
        return n*conv_ato_cm2*2*np.sum(cache)/len(cache)
    l=0
    for n in [2,3,4,5]:
        photo = []
        RN    = find_RN(n)
        Epp   = np.linspace(conv*0.51/n**2,conv*10/27.2114,51)
        for Ephoto in Epp:
            photo.append(cross_section_photoionization(l,n,Ephoto))
        photoionization = CubicSpline(Epp,photo)
        def fc_photoionization(E,n):
            if E < conv*0.51/n**2:
                return 0
            else:
                return photoionization(E)
        
        energies = np.logspace(np.log(Ec[n]),np.log(10),31,base=np.exp(1))
        RE_list  = [find_RE(E) for E in energies]
        sigma    = Parallel(n_jobs=num_cores)(delayed(fc_sigma)(energies[i],n,RN,RE_list[i]) for i in trange(len(energies)))
        #plt.plot(energies,sigma,color=color[n],ls='-',label='$n='+str(n)+'$',zorder=11)
        #print(n)
        #print(sigma)
        np.savetxt('output/alt_sigma_n='+str(n)+'_l='+str(l)+'.txt',np.transpose([energies,sigma]))
     
    '''
    array = pd.read_csv('H(1s)-H(2s).csv')
    xdata = [float(i) for i in array['xdata']]
    ydata = [float(i) for i in array['ydata']]
    plt.plot(xdata,ydata,'-*',color='yellow',label='H(1s)+H(2s)')

    array = pd.read_csv('H(1s)-H(3s).csv')
    xdata = [float(i) for i in array['xdata']]
    ydata = [float(i) for i in array['ydata']]
    plt.plot(xdata,ydata,'-x',color='black',label='H(1s)+H(3s)')

    array = pd.read_csv('H(1s)-H(4s).csv')
    xdata = [float(i) for i in array['xdata']]
    ydata = [float(i) for i in array['ydata']]
    plt.plot(xdata,ydata,'-+',color='olive',label='H(1s)+H(4s)')
    '''

    '''
    array = pd.read_csv('n_10.csv')
    xdata = [float(i)**2 for i in array['xdata']]
    ydata = [float(i) for i in array['ydata']]
    plt.plot(xdata,ydata,'-',color='grey')

    array = pd.read_csv('n_15.csv')
    xdata = [float(i)**2 for i in array['xdata']]
    ydata = [float(i) for i in array['ydata']]
    plt.plot(xdata,ydata,'-',color='grey')

    array = pd.read_csv('n_20.csv')
    xdata = [float(i)**2 for i in array['xdata']]
    ydata = [float(i) for i in array['ydata']]
    plt.plot(xdata,ydata,'-',color='grey')

    array = pd.read_csv('n_25.csv')
    xdata = [float(i)**2 for i in array['xdata']]
    ydata = [float(i) for i in array['ydata']]
    plt.plot(xdata,ydata,'-',color='grey')
    '''

    '''
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$E$ [eV]')
    plt.ylabel('$\\sigma(E)$ [$10^{-16}$ cm$^2$]')
    plt.xlim(1e-2,1e1)
    plt.legend()
    plt.savefig('sigma.png',dpi=300)
    '''
    
    










    print("--- %s seconds ---" % (time.time() - start_time))
