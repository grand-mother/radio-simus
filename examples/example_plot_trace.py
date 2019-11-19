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


TODO: add a SNR calculation as example how to treat signal
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
from astropy import units as u  
import h5py

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
import radio_simus 
from radio_simus.utils import _getAngle
from radio_simus.signal_treatment import p2p, hilbert_env, hilbert_peak
from radio_simus.io_utils import _load_to_array, _load_eventinfo_fromhdf


if __name__ == '__main__':

    if ( len(sys.argv)<2 ):
        print("""
        Script to read in a single antenna trace
        -- Plot trace (PLOT)
        -- Plot Hilbert envelope (HILB)
        -- calculate Angle (ANGLE) -- to be fixed for CR (modules.py)
        
        Usage: python3 example_plot_trace.py <path to event> <keyword:efield/voltages/...> <antID>
        example: python3 example_plot_trace.py event_000001.hdf5 efield 9

        """)
        sys.exit(0)

# path to event/hdf5 file
path_hdf5 = sys.argv[1]
# antenna ID
i = sys.argv[3]

keyword=str(sys.argv[2]) # efield / voltages

# load trace
trace, time_unit, unit, ant_position, ant_slopes = _load_to_array(path_hdf5, content=keyword, ant="/" + str(i) +"/" )
time_unit=str(time_unit)
unit=str(unit)





########### How to plot a trace 
PLOT=True
if PLOT:
    #### Plotting section 
    fig=plt.figure(1,figsize=(8, 5), dpi=120, facecolor='w', edgecolor='k')

    plt.title(keyword +", antenna ID: "+str(i))
    plt.plot(trace.T[0], trace.T[1], 'b-', label= "NS, PtP = {0:.2f} ".format(p2p(trace.T[1]))+unit )
    plt.plot(trace.T[0], trace.T[2], 'r-', label= "EW, PtP = {0:.2f} ".format(p2p(trace.T[2]))+unit )
    plt.plot(trace.T[0], trace.T[3], 'g-', label= "Up, PtP = {0:.2f} ".format(p2p(trace.T[3]))+unit )

    plt.xlabel(r"time ("+time_unit+")", fontsize=16)
    plt.ylabel(r"Amplitude ("+unit+")", fontsize=16)
    plt.legend(loc='best', fancybox=True)

    plt.show()
    
########### Hilbert envelope, maximum and its time
HILB=False
if HILB:
    hilbert = hilbert_env(trace.T[2])
    print("peak amplitude and its time :",hilbert_peak(trace.T[0], trace.T[2]))
    
    fig2=plt.figure(2,figsize=(8, 5), dpi=120, facecolor='w', edgecolor='k')

    plt.title("Hilbert env -- antenna ID: "+str(i))
    plt.plot(trace.T[0], hilbert, 'b-', label= "EW, Max.= {0:.2f} ".format(max(hilbert))+unit )
    plt.plot(trace.T[0], trace.T[2], 'r--', label= "EW, PtP = {0:.2f} ".format(p2p(trace.T[2]))+unit )

    plt.xlabel(r"time ("+time_unit+")", fontsize=16)
    plt.ylabel(r"Amplitude ("+unit+")", fontsize=16)
    plt.legend(loc='best', fancybox=True)

    plt.show()

########### getAngle between shower axis and antenna position -- only works for neutrino induced showers for now
ANGLE=False
if ANGLE: # Only for neutrino showers
    from radio_simus.modules import _getXmax, _dist_decay_Xmax, _get_XmaxPosition, _get_CRzenith
    
    #load shower information
    shower, antID, positions, slopes = _load_eventinfo_fromhdf(path_hdf5)
    
    # this does not yet work correctly....
    Xmaxpos =  _get_XmaxPosition(shower["primary"], shower["energy"].value, 
                                 shower["zenith"].value, shower["azimuth"].value, shower["injection_height"].value)
    print(shower, Xmaxpos)
    #zen_CR = _get_CRzenith(shower["zenith"].value,shower["injection_height"].value, 0.)
    print("Angle in deg: ", _getAngle(refpos=Xmaxpos, theta=shower["zenith"].value, azim=shower["azimuth"].value, 
                                      ANTENNAS=positions[np.where(antID==str(i))].value, core=shower["core"].value))
    
    ## explicit example
    #Xmaxpos =  _get_XmaxPosition('electron', 0.5*1e18, 87., 0., 2800.)
    #print(_getAngle(refpos=Xmaxpos ,theta=87.,azim=0,ANTENNAS=np.array([50051,0., 3272]), core=[67407.46,0.,6332.68]))
    


