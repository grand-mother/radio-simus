'''
Script to read in a whole event in hdf5 files and to plot the field strengths.


Usage: python3 example_plot_2D.py <path to event in hdf5 format> <input format>
Example: python3.6 example_plot_2D.py event_000001.hdf5 v

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
    
TODO: How to handle not-existing positions    
    
'''

import sys
from sys import argv
import os
import glob

import numpy as np
from numpy import *

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
#from mpl_toolkits.mplot3d import Axes3D ## leading to troubles with python3

#from astropy import units as u
#muV_m = u.u * u.V / u.m

from astropy.table import Table
from astropy.table import hstack
import h5py

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))

# load config file first
import radio_simus 
radio_simus.load_config('./test.config')
arrayfile = radio_simus.config.arrayfile

from radio_simus.signal_treatment import p2p
from radio_simus.io_utils import _load_path,_load_to_array, _load_eventinfo_fromhdf



### path to full antenna array to plot in the background
ant_array = np.loadtxt(arrayfile,comments="#")
print("Loading full array: ", arrayfile)


def plot_2D(path1, key, save=None):
    ''' Creates a 2D plot of the signal distribution in the array from hdf5 event file
    '''
    

    ## get information on event from hdf5 file -- in principle not needed here
    shower, ant_ID,  positions, slopes = _load_eventinfo_fromhdf(path1)
    #print( len(positions[0]) , "antenna positions simulated in Event ", shower["ID"])
    
        
    ## get information on analysis already performed
    analysis, ana_info, ana_meta = _load_path(path1, path="/analysis")
    
    # short way -- but works only for voltages since the p2p values in hdf5 files are only valid for voltage traces by default
    if key == "voltages":
        trigger_any=analysis["trigger_aggr_any"]
        trigger_xy=analysis["trigger_aggr_xy"]
        p2p_Ex=analysis["p2p_x"]
        p2p_Ey=analysis["p2p_y"]
        p2p_Ez=analysis["p2p_z"]
        p2p_total=analysis["p2p_xy"]
        unit=str(analysis["p2p_x"].unit)
        
        x_pos=positions.T[0]
        y_pos=positions.T[1]
        z_pos=positions.T[2]
        
    # long way
    if key == "efield":   
        # create an array of simulated observer positions
        p2p_Ex = []
        p2p_Ey = []
        p2p_Ez = []
        p2p_total = []

        x_pos=[]
        y_pos=[]
        z_pos=[]

        trigger_any=[]
        trigger_xy=[]
        

        f = h5py.File(path1, 'r') # open event files
        for ID in f.keys(): #loop over antennas in event 
            try:
                
                # Note: Be aware of that trigger was done on voltages
                if analysis["trigger_aggr_xy"][np.where(ant_ID==ID)]:
                    trigger_xy.append(ID)
                if analysis["trigger_aggr_any"][np.where(ant_ID==ID)]:
                    trigger_any.append(ID)
                    
                # obtain electric field trace from hdf5 file
                trace, time_unit1, unit1, ant_position, ant_slopes = _load_to_array(path1, content=key, ant=ID)
                unit=str(unit1)
                time_unit=str(time_unit1)  
                
                # get antenna position
                #position = positions[np.where(ant_ID==ID)].value[0]
                position = ant_position
                x_pos.append(position[0].value)
                y_pos.append(position[1].value)
                z_pos.append(position[2].value)   
                         
                ## call p2p function 
                ### or mcuh simpler: call p2p from analysis in hdf outside the loop
                p2ps = p2p(trace)
                
                # now it depends of which parameter we need the peak amplitude: here we go for the peak-to-peak amplitude
                p2p_Ex.append(p2ps[0])
                p2p_Ey.append(p2ps[1])
                p2p_Ez.append(p2ps[2]) 
                p2p_total.append(p2ps[4])

            except: # in case file does not exist
                continue
                

    #################
    #PLOTTING SECTION
    #################

    #c_array = np.zeros(len(ant_array.T[0]))

    ###### Plot a 3d figure of the total peak amplitude to see the actual array and the the signal distribution ## leading to troubles with python3
    #fig = plt.figure(1, figsize=(8,7), dpi=120, facecolor='w', edgecolor='k')
    #ax = fig.add_subplot(111, projection='3d')

    #ax.scatter(ant_array.T[1], ant_array.T[2],ant_array.T[3] , c="None",  
            #marker='o', edgecolor='gray', cmap=cm.gnuplot2_r)
    #col=ax.scatter(x_pos, y_pos,z_pos , c=p2p_total,  vmin=min(p2p_total), vmax=max(p2p_total),  
                #marker='o', edgecolor='black', cmap=cm.gnuplot2_r)
    #cbar=plt.colorbar(col)
    #cbar.set_label('peak-to-peak amplitude distribution in '+unit, labelpad=+1)
    #try:
        #core=txt.meta["core"]
        #ax.scatter(core[0], core[1],core[2] ,c="red", marker='x', s=50)
    #except:
        #print("no core info available")

    #ax.set_xlabel('positions along NS (m)')
    #ax.set_ylabel('positions along EW (m)')
    #ax.set_zlabel('positions Up (m)')
    ##plt.title('peak-to-peak amplitude distribution in '+unit)



    ##### Plot a 2d figures of the NS, ES, UP component and total peak amplitude in positions along North-South and East-West 
    fig2 = plt.figure(2,figsize=(9,7), dpi=120, facecolor='w', edgecolor='k')
    
    ax1=fig2.add_subplot(221)
    plt.title("Event: "+ str(shower["ID"])+"\n #triggered: any ="+ str(len(trigger_any)) + " xy = "+ str(len(trigger_xy)), fontsize=10)
    name = 'NS-component ('+unit+')'
    #plt.title(name)
    ax1.set_xlabel('positions along NS (m)')
    ax1.set_ylabel('positions along EW (m)')
    ax1.scatter(ant_array.T[1], ant_array.T[2], c="None",  
            marker='o', edgecolor='gray', cmap=cm.gnuplot2_r)
    col1=ax1.scatter(x_pos, y_pos, c=p2p_Ex,  vmin=min(p2p_Ex), vmax=max(p2p_Ex),  marker='o', edgecolor='black', cmap=cm.gnuplot2_r)
    cbar1=plt.colorbar(col1)
    cbar1.set_label(name, labelpad=+1)
    try:
        core=txt.meta["core"]
        ax1.scatter(core[0], core[1] ,c="red", marker='x', s=50)
    except:
        print("no core info available")

    plt.tight_layout()
        
    ax2=fig2.add_subplot(222)
    name = 'EW-component ('+unit+')' 
    #plt.title(name)
    ax2.set_xlabel('positions along NS (m)')
    ax2.set_ylabel('positions along EW (m)')
    ax2.scatter(ant_array.T[1], ant_array.T[2], c="None",  
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
    ax3.scatter(ant_array.T[1], ant_array.T[2], c="None",  
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
    ax4.scatter(ant_array.T[1], ant_array.T[2], c="None",  
            marker='o', edgecolor='gray', cmap=cm.gnuplot2_r)
    col4=ax4.scatter(x_pos, y_pos, c=p2p_total,  vmin=min(p2p_total), vmax=max(p2p_total),  marker='o', edgecolor='black', cmap=cm.gnuplot2_r)
    cbar1=plt.colorbar(col4)
    cbar1.set_label(name, labelpad=+1)
    try:
        ax4.scatter(core[0], core[1] ,c="red", marker='x', s=50)
    except:
        print("no core info available")
        
    plt.tight_layout()

    if save!=None:
        plt.savefig(save)
    else:    
        plt.show()


def main():
    if ( len(sys.argv)<3 ):
        print("""
            Example on how to do a 2D plot, show the radio footprint in the array
                -- read in original array list (from config-file)
                -- plots the p2p distribution componentwise
                -- reads in single antenna hdf5 files
            
            Usage: python3 example_plot_2D.py <path to event folder> <e/v>(efield or voltages to plot)
            Example: python3 example_plot_2d.py ../../CoREAS/GP300_test2/000001/ v
            
        """)
        sys.exit(0)
        
  # path to folder containing the inp-file, trace files and antpos.dat 
    path = sys.argv[1]

    if sys.argv[2] == 'e' or sys.argv[2] == 'efield':
        key = "efield" 
    if sys.argv[2] == 'v' or sys.argv[2] == 'voltages':
        key = "voltages" 

        
        
    plot_2D(path, key, save=None)
  
if __name__== "__main__":
  main()


    
