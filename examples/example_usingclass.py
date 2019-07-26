
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
    index = np.where(positions == pos[:])[0]
    #print('indexx', index, positions[index], pos)
    return antIDs[index]

################################


# path to folder containing the inp-file, trace files and antpos.dat 
path = sys.argv[1]


# create an array
p2p_Ex = []
p2p_Ey = []
p2p_Ez = []
p2p_total = []

x_pos=[]
y_pos=[]
z_pos=[]

i=0

#create  "empty detector"
det = detector()
#create detector from file defined in config file
create_from_file(det, arrayfile)

# get all antennas positions
array = get_array(det)
# get all slopes
slopes = get_slopes(det)

antIDs = array.T[0]
positions = array[:,1:] 


# loop over all folder


for file in glob.glob(path+"*.hdf5"):
    f = Table.read(file, path='efield') 
    print("\n simulated position ", f.meta["position"])
    
    ## find antenna position and its slope per ID - works
    #ID = int(file.split('/')[-1].split('.hdf5')[0].split('table_')[-1])
    #pos_ant = find_antennaposition(det, ID)
    #pos_slope = find_antennaslope(det, ID)
    #print("----", ID, pos_ant, pos_slope)
    
    
    pos_ID = find_antennaID(positions, antIDs, f.meta["position"])
    ##print("----", f.meta["position"], pos_ID)
    
    if i==0: # just get the first antenna to readin meta info
        testshower = sim_shower()
        loadInfo_toShower(testshower, f.meta)
        #print(testshower.get_all())
    i+=1
    
    try:
        g = Table.read(file, path='voltages') 
    except IOError:
        print("voltages not computed")
    
    # check whether voltages exists, if not compute voltage
    
    
    # check whether antenna ID position and slope already exits, otherwise load to detector
    


# plot full array (gray), simulated position with voltages and mark triggers (red) for each event, save as png
