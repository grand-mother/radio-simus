'''
Scripts reads in one hdf5file == one antenna and plots it


python3 example_plot_2D.py <path> <input format> <antenna ID> 
Parameters:
path: path to folder containing the trace files and antpos.dat
input format:   e: readin electric field traces, 
                #####f: read in filtered traces (f1, f2 in MHz hardcoded in the script, 
                v: read in voltage traces 
antenna ID: number of desired antenna
                
Output:
will plot all 3 field components for an antenna 
Note: This is just an example file for reading-in and plotting signals of a full array.
It is far from being perfect, just for beginners in the hand-on session
'''

### frequency, used if 'f' is chosen in MHz
f1 = 50
f2 = 200


import sys
from sys import argv
import os

import numpy as np
from numpy import *

import matplotlib.pyplot as plt
import pylab

from astropy.table import Table
import h5py

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
import radio_simus 
from radio_simus.utils import p2p


# path to folder containing the inp-file, trace files and antpos.dat 
path = sys.argv[1]
# antenna ID
i = sys.argv[3]

# load trace
if sys.argv[2] == 'e': # readin electric field trace
    txt = Table.read(path+ 'table_a'+str(i)+'.hdf5', path="efield") #txt.T[0]:time in ns, txt.T[1]: North-South, txt.T[2]: East-West, txt.T[3]: Up , all electric field in muV/m
    trace = np.array([txt['Time'], txt['Ex'], txt['Ey'], txt['Ez']]).T
    print(trace)
    unit='muV/m'
if sys.argv[2] == 'v': # readin voltage trace
    txt = Table.read(path+ 'table_a'+str(i)+'.hdf5', path="voltages")#txt.T[0]: time in ns, txt.T[1]: North-South, txt.T[2]: East-West, txt.T[3]: Up , all voltage in muV
    trace=np.array([txt['Time'], txt['Vx'], txt['Vy'], txt['Vz']]).T
    unit='muV (' +str(f1)+'-'+str(f2)+'MHz)'
#if sys.argv[2] == 'f': # readin filtered electric field trace
    #txt = np.loadtxt(path+ 'a'+str(i)+'_'+str(f1)+'-'+str(f2)+'MHz.dat') #txt.T[0]: time in ns, txt.T[1]: North-South, txt.T[2]: East-West, txt.T[3]: Up , all electric field in muV/m
    #unit='muV/m (' +str(f1)+'-'+str(f2)+'MHz)'



#### Plotting section 
fig=plt.figure(1,figsize=(8, 5), dpi=120, facecolor='w', edgecolor='k')

plt.title("antenna ID: "+str(i))
plt.plot(trace.T[0], trace.T[1], 'b-', label= "Ex=NS, PtP={0:.2f}".format(p2p(trace.T[1]))+unit )
plt.plot(trace.T[0], trace.T[2], 'r-', label= "Ey=EW, PtP={0:.2f}".format(p2p(trace.T[2]))+unit )
plt.plot(trace.T[0], trace.T[3], 'g-', label= "Ez=Up, PtP={0:.2f}".format(p2p(trace.T[3]))+unit )

plt.xlabel(r"time (ns)", fontsize=16)
plt.ylabel(r"Amplitude ("+unit+")", fontsize=16)
plt.legend(loc='best', fancybox=True)

plt.show()
