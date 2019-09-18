'''
    Example on how to use classes
    -- Analysis trigger for events, create a list of events (class objects) and adds trigger info 
    
    Use: python3 example_usingclass.py <folder event set>
    Example: python3 example_usingclass.py ../../CoREAS/GP300_centered/
    
    NOTE: still ongoing, already usable

'''



import sys
from sys import argv
import os
import glob

import numpy as np
from numpy import *

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
#from mpl_toolkits.mplot3d import Axes3D

#from astropy import units as u
#muV_m = u.u * u.V / u.m

from astropy.table import Table
from astropy.table import hstack
import h5py

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
#import radio_simus 
from radio_simus.signal_treatment import p2p
from radio_simus.shower import *
from radio_simus.detector import detector, create_from_file, get_array, get_slopes, find_antennaposition, find_antennaslope
from radio_simus.__init__ import arrayfile 



##############################
### define new functions and test

#def find_antennaposition(positions, antIDs, ID):
    #index = np.where(antIDs == int(ID))[0][0]
    ##print(index)
    #return positions[index]

def find_antennaID(positions, antIDs, pos):
    pos=np.around(pos)
    print(type(positions), type(pos))
    print(positions[0], pos)
    index = np.where((positions[:]== pos).all())
    print('indexx', index, positions[index,0])
    return antIDs[index]
## not yet working
#pos_ID = find_antennaID(positions, antIDs, np.asarray(f.meta["position"]))
#print("----", f.meta["position"], pos_ID, "original", ID)
################################


# path to folder containing the inp-file, trace files and antpos.dat 
eventfolder = sys.argv[1]


# create an array
p2p_Ex = []
p2p_Ey = []
p2p_Ez = []
p2p_total = []

x_pos=[]
y_pos=[]
z_pos=[]


#create  "empty detector"
det = detector()
#create detector=antenna array from file defined in config file
create_from_file(det, arrayfile)

# get all antennas positions
array = get_array(det)
# get all slopes
slopes = get_slopes(det)

antIDs = array.T[0]
positions = array[:,1:4] 

event = []

# loop over all folder
for path in glob.glob(eventfolder+"/*"):
    if os.path.isdir(path): # only pick event folders
        print("\n==> Reading Event from:", path)

    
        # loop over all antenna positions in event
        i=0
        trigger_any=[]
        trigger_xy=[]
        

        for file in glob.glob(path+"/*.hdf5"):
            f = Table.read(file, path='efield') 
            #print("\n simulated position ", f.meta["position"])
            
            ## find antenna position and its slope per ID - works
            ID = int(file.split('/')[-1].split('.hdf5')[0].split('table_')[-1])
            #pos_ant = find_antennaposition(det, ID)
            #pos_slope = find_antennaslope(det, ID)

            
            if i==0: # just get the first antenna to readin meta info
                testshower = sim_shower()
                loadInfo_toShower(testshower, f.meta)
                param = testshower.get_all() # get all parameters, all call them separately
                print("=== SUMMARY EVENT ===")
                print("ShowerID = ",  param[0], " primary = ", param[1], " energy/eV = ", param[2] , " zenith/deg = ", param[3], " azimuth/deg = ", param[4],  " injectionheight/m = ", param[5] )
                
                event.append(testshower)
            i+=1
            
            
            ## read voltages for analysis
            try:
                g = Table.read(file, path='voltages') 
                if g.meta["trigger"][0] ==1:
                    trigger_any.append(ID)
                if g.meta["trigger"][1] ==1:
                    trigger_xy.append(ID)
            except IOError:
                print("voltages not computed for antenna: ", str(ID))
            
            # check whether voltages exists, if not compute voltage
            
            
            # check whether antenna ID position and slope already exits, otherwise load to detector
            
        ## Trigger Analysis
        if len(trigger_any)>5 or len(trigger_xy)>5:
            print(" ===> shower would have triggered: any =", len(trigger_any), " xy = ", len(trigger_xy))
            # TODO: add trigger info to class
            print(event[-1].showerID)
            event[-1].add_trigger(1)
            
        else:
            testshower.add_trigger(0)
        
            
        # plot full array (gray), simulated position with voltages and mark triggers (red) for each event, save as png
        ## could be similar to example_plot_2D
    
    else: 
        continue
    
print("How to handle now the list of events...")    
print(event[0].showerID, event[0].trigger)
