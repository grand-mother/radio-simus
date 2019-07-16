'''
Script to read in a whole event in hdf5 files and to plot the field strengths.


Usage: python3 example_plot_2D.py <path to event in hdf5 format> <input format>
Example: python3.6 example_plot_2D.py ../../CoREAS/GP300-v3/000002/ v

Parameters:
    path: path to folder containing the trace files and antpos.dat
    input format:   e: readin electric field traces, 
                    ####f: read in filtered traces (f1, f2 in MHz hardcoded in the scripr, 
                    v: read in voltage traces 
Output:
    - will plot the total peak-to-peak distribution for an array in 3D
    - will plot the peak-to-peak distribution for the single components as a 2D scatter plot
Note: This is just an example file for reading-in and plotting signals of a full array.
    It is far from being perfect, just for beginners in the hand-on session
'''

import sys
from sys import argv
import os
import glob

import numpy as np
from numpy import *

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
from mpl_toolkits.mplot3d import Axes3D

#from astropy import units as u
#muV_m = u.u * u.V / u.m

from astropy.table import Table
from astropy.table import hstack
import h5py

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
import radio_simus 
from radio_simus.utils import p2p



# path to folder containing the inp-file, trace files and antpos.dat 
path = sys.argv[1]

if sys.argv[2] == 'e' or sys.argv[2] == 'efield':
    #txt.T[0]: time in ns, txt.T[1]: North-South, txt.T[2]: East-West, txt.T[3]: Up , all electric field in muV/m
    key = "efield" 
    unit=r'$\mu$V/m'
#if sys.argv[2] == 'f':
    ##txt.T[0]: time in ns, txt.T[1]: North-South, txt.T[2]: East-West, txt.T[3]: Up , all voltage in muV
    #ending = 'out_*.txt'
    #unit=r'$\mu$V/m (' + str(f1)+'-'+str(f2)+'MHz)'
if sys.argv[2] == 'v' or sys.argv[2] == 'voltages':
    #txt.T[0]: time in ns, txt.T[1]: North-South, txt.T[2]: East-West, txt.T[3]: Up , all electric field in muV/m
    key = "voltages" 
    unit=r'$\mu$V'
    

### path to full antenna array to plot in the background
#ant_array = np.loadtxt('/home/laval1NS/zilles/CoREAS/AntennasAnneG300.txt',comments="#")
#print("Loading full array: ", '/home/laval1NS/zilles/CoREAS/AntennasAnneG300.txt')
ant_array = np.loadtxt('/home/laval1NS/zilles/CoREAS/regular_array.txt',comments="#")
print("Loading full array: ", '/home/laval1NS/zilles/CoREAS/regular_array.txt')

# create an array of simulated observer positions
p2p_Ex = []
p2p_Ey = []
p2p_Ez = []
p2p_total = []

x_pos=[]
y_pos=[]
z_pos=[]

for file in glob.glob(path+"*.hdf5"):   
    # get ID and antenna position
    base=os.path.basename(file)
    # coreas
    ID=(os.path.splitext(base)[0]).split("_a")[1] # remove table_a 

    try: # load corresponding file
        txt = Table.read(file, path=key)
        if key == "efield":
            keywords = [txt['Time'], txt['Ex'], txt['Ey'], txt['Ez']]
            ## TODO check units - no idea zet to hanlde it properly
            print(" Is efield in muV/m: ",txt['Ex'].unit)
        if key == "voltages":
            keywords = [txt['Time'], txt['Vx'], txt['Vy'], txt['Vz']]
            ## TODO check units - no idea zet to hanlde it properly
            print(" Is voltages in muV: ",txt['Vx'].unit)
        
        
        # get antenna position
        position = txt.meta['position']
        x_pos.append(position[0])
        y_pos.append(position[1])
        z_pos.append(position[2])

                
        # convert astropy table to np array
        trace=np.array(keywords).T

        ### no solution find yet to treat the units nicely
        #txt.T[0]*=u.nanosecond 
        #txt.T[1]*=muV_m
        #txt.T[2]*=muV_m
        #txt.T[3]*=muV_m

        ### call p2p function
        p2ps = p2p(trace)
        # now it depends of which parameter we need the peak amplitude: here we go for the peak-to-peak amplitude
        p2p_Ex.append(p2ps[0])
        p2p_Ey.append(p2ps[1])
        p2p_Ez.append(p2ps[2]) 
        
        amplitude = (np.sqrt(trace.T[1]**2. + trace.T[2]**2. + trace.T[3]**2.)) # combined components
        p2p_total.append(p2p(amplitude))

    except IOError: # in case file does not exist
        print("-- File ", base, " does not exist -> set to 0.")
        p2p_Ex.append(0.)
        p2p_Ey.append(0.)
        p2p_Ez.append(0.)
        p2p_total.append(0.)
        
    #print(x_pos[-1].to(u.m))




#################
#PLOTTING SECTION
#################

c_array = np.zeros(len(ant_array.T[0]))

##### Plot a 3d figure of the total peak amplitude to see the actual array and the the signal distribution
fig = plt.figure(1, figsize=(8,7), dpi=120, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection='3d')

ax.scatter(ant_array.T[0], ant_array.T[1],ant_array.T[2] , c="None",  
           marker='o', edgecolor='gray', cmap=cm.gnuplot2_r)
col=ax.scatter(x_pos, y_pos,z_pos , c=p2p_total,  vmin=min(p2p_total), vmax=max(p2p_total),  
               marker='o', edgecolor='black', cmap=cm.gnuplot2_r)
cbar=plt.colorbar(col)
cbar.set_label('peak-to-peak amplitude distribution in '+unit, labelpad=+1)
try:
    core=txt.meta["core"]
    ax.scatter(core[0], core[1],core[2] ,c="red", marker='x', s=50)
except:
    print("no core info available")

ax.set_xlabel('positions along NS (m)')
ax.set_ylabel('positions along EW (m)')
ax.set_zlabel('positions Up (m)')
#plt.title('peak-to-peak amplitude distribution in '+unit)



##### Plot a 2d figures of the NS, ES, UP component and total peak amplitude in positions along North-South and East-West 
fig2 = plt.figure(2,figsize=(12,7), dpi=120, facecolor='w', edgecolor='k')
    
ax1=fig2.add_subplot(221)
name = 'NS-component ('+unit+')'
#plt.title(name)
ax1.set_xlabel('positions along NS (m)')
ax1.set_ylabel('positions along EW (m)')
ax1.scatter(ant_array.T[0], ant_array.T[1], c="None",  
           marker='o', edgecolor='gray', cmap=cm.gnuplot2_r)
col1=ax1.scatter(x_pos, y_pos, c=p2p_Ex,  vmin=min(p2p_Ex), vmax=max(p2p_Ex),  marker='o', edgecolor='black', cmap=cm.gnuplot2_r)
cbar1=plt.colorbar(col1)
cbar1.set_label(name, labelpad=+1)
try:
    ax1.scatter(core[0], core[1] ,c="red", marker='x', s=50)
except:
    print("no core info available")

plt.tight_layout()
    
ax2=fig2.add_subplot(222)
name = 'EW-component ('+unit+')' 
#plt.title(name)
ax2.set_xlabel('positions along NS (m)')
ax2.set_ylabel('positions along EW (m)')
ax2.scatter(ant_array.T[0], ant_array.T[1], c="None",  
           marker='o', edgecolor='gray', cmap=cm.gnuplot2_r)
col2=ax2.scatter(x_pos, y_pos, c=p2p_Ey,  vmin=min(p2p_Ey), vmax=max(p2p_Ey),  marker='o', edgecolor='black', cmap=cm.gnuplot2_r)
cbar1=plt.colorbar(col2)
cbar1.set_label(name, labelpad=+1)
try:
    ax2.scatter(core[0], core[1] ,c="red", marker='x', s=50)
except:
    print("no core info available")
    
plt.tight_layout()
    
ax3=fig2.add_subplot(223)
name = 'Up-component ('+unit+')' 
#plt.title(name)
ax3.set_xlabel('positions along NS (m)')
ax3.set_ylabel('positions along EW (m)')
ax3.scatter(ant_array.T[0], ant_array.T[1], c="None",  
           marker='o', edgecolor='gray', cmap=cm.gnuplot2_r)
col3=ax3.scatter(x_pos, y_pos, c=p2p_Ez,  vmin=min(p2p_Ez), vmax=max(p2p_Ez),  marker='o', edgecolor='black', cmap=cm.gnuplot2_r)
cbar1=plt.colorbar(col3)
cbar1.set_label(name, labelpad=+1)
try:
    ax3.scatter(core[0], core[1] ,c="red", marker='x', s=50)
except:
    print("no core info available")
    
plt.tight_layout()

ax4=fig2.add_subplot(224)
name = 'total ('+unit+')'
#plt.title(name)
ax4.set_xlabel('positions along NS (m)')
ax4.set_ylabel('positions along EW (m)')
ax4.scatter(ant_array.T[0], ant_array.T[1], c="None",  
           marker='o', edgecolor='gray', cmap=cm.gnuplot2_r)
col4=ax4.scatter(x_pos, y_pos, c=p2p_total,  vmin=min(p2p_total), vmax=max(p2p_total),  marker='o', edgecolor='black', cmap=cm.gnuplot2_r)
cbar1=plt.colorbar(col4)
cbar1.set_label(name, labelpad=+1)
try:
    ax4.scatter(core[0], core[1] ,c="red", marker='x', s=50)
except:
    print("no core info available")
    
plt.tight_layout()


plt.show()
