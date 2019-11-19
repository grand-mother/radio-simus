################################
#### by A. Zilles
################################
#!/usr/bin/env python

import sys

from sys import argv
import os
import glob

import time

import numpy as np
from numpy import *


import logging
logging.basicConfig(filename="example_usingclass.log", level=logging.INFO)
print("Log-file produced: ", "example_usingclass.log")
logger = logging.getLogger('Main')

import tqdm

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
#from radio_simus.detector import detector, create_from_file, get_array, get_slopes, find_antennaposition, find_antennaslope
from radio_simus.detector import Detector
from radio_simus.__init__ import arrayfile 




if __name__ == '__main__':
    
    if ( len(sys.argv)<2 ):
        print("""
        Example on how to use classes
        -- Analysis trigger for events, create a list of events (class objects) and  trigger 1/0 to class attributes
        -- create a png with statistic for triggering
        
        Use: python3 example_usingclass.py <folder event set>
        Example: python3 example_usingclass.py ../../CoREAS/GP300_centered/
    
        """)
        sys.exit(0)
    
    
    

    # path to folder containing the inp-file, trace files and antpos.dat 
    eventfolder = sys.argv[1]


    ### FANCY THINGS YOU CAN DO WITH THE DETECTOR CLASS
    ### but not needed for this short study
    ##create  "empty detector"
    det = Detector()
    #create detector=antenna array from file defined in config file
    det.create_from_file(arrayfile)

    ## get all antennas positions
    #array = det.position
    ## get all slopes
    #slopes = det.slope
    ## antenna IDs
    #antIDs = det.ID

    #print(det.position)
    #print(det._attributes, det.origin, det.location)
    ##print(det.position.T, det.ID)
    #print(det.find_position(1), det.find_slope(1))
    #print(det.find_antenna(1))
            
    #ant=det.find_antenna(1)
    #print(ant.position)




    print("\nScan of events ...")
    
    event = [] # python list        
        
    # loop over all folder
    #for path in tqdm.tqdm(glob.glob(eventfolder+"/*/")):
    for path in glob.glob(eventfolder+"*.hdf5"): # loop over files in folder
        if os.path.isfile(path): # only pick event folders
            logger.debug("... Reading Event from:"+ path)


            # create shower object and set attributes
            testshower = SimulatedShower()
            g=Table.read(path, path="/event")
            loadInfo_toShower(testshower, g.meta)
            #shower, positions, slopes = _load_eventinfo(file)
            
            logger.info("   SUMMARY EVENT: ShowerID = "+  str(testshower.showerID)
                        + " primary = "+ str(testshower.primary)+ " energy/eV = "+ str(testshower.energy) 
                        + " zenith/deg = "+ str(testshower.zenith)+ " azimuth/deg = "+ str(testshower.azimuth)
                        + " injectionheight/m = "+ str(testshower.injectionheight) )
                        
            event.append(testshower)

            
            analysis=Table.read(path, path="/analysis")
            #print(analysis["trigger_aggr_xy"], analysis["trigger_aggr_any"], analysis["trigger_cons_xy"],analysis["trigger_cons_any"])
             
            ## EXAMPLE: Trigger Analysis
            if np.sum(analysis["trigger_aggr_xy"])>5 or np.sum(analysis["trigger_aggr_any"])>5:
                logger.info("   => shower triggers (aggr): any =" + str(np.sum(analysis["trigger_aggr_xy"])) 
                            + " xy = " + str(np.sum(analysis["trigger_aggr_any"])))
                # add trigger info to class
                event[-1].trigger=1
            else:
                event[-1].trigger=0
                
                    
            ## TODO plot full array (gray), simulated position with voltages and mark triggers (red) for each event, save as png
            ### could be similar to example_plot_2D
            
        #else: 
            #continue


    ###### START ANALYSIS ################### 
    print("\nStart an analysis ...")    

    # print attributes - define by object class, does not mean that they are !=None
    print("\nAvailable attributes: ", testshower._attributes ) # refers to last event
    # --> Available attributes:  ('showerID', 'primary', 'energy', 'zenith', 'azimuth', 'injectionheight', 'trigger', 'simulation', 'Xmax')

    Event_ID=list(map(lambda i: i.showerID, event))

    ### Calculate Ratio of detected events
    trigger=list(map(lambda i: i.trigger, event))
    print("\n"+str(sum(trigger))+" out of "+str(len(trigger))+" events detected --> "+str(100.* sum(trigger)/len(trigger))+"% detection rate"+"\n")    

    ### find triggered events
    trigger=np.asarray(trigger)
    index = np.where(trigger==1)[0]

    # parameters
    energy=np.asarray(list(map(lambda i: i.energy/u.eV, event)))
    zenith=np.asarray(list(map(lambda i: i.zenith/u.deg, event)))
    azimuth=np.asarray(list(map(lambda i: i.azimuth/u.deg, event)))
    primary=np.asarray(list(map(lambda i: i.primary, event)))





    ##### PLOTTING
    # Plot
    plt.rcParams.update({'figure.figsize':(12,5)})

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
    ax1.hist(energy, color="orange", bins=50, alpha=0.5, label='all = '+str(len(trigger)))
    ax1.hist(energy[index], bins=50, label='trigger (aggr.) = '+str(sum(trigger)))
    ax1.legend()
    ax1.set_title('energy in eV')
    ax2.hist(zenith, color="orange", bins=50, alpha=0.5)
    ax2.hist(zenith[index], bins=50)
    ax2.set_title('zenith in deg')
    ax3.hist(azimuth, color="orange", bins=50, alpha=0.5)
    ax3.hist(azimuth[index], bins=50)
    ax3.set_title('azimuth in deg')
    ax4.hist(primary, color="orange", bins=50, alpha=0.5)
    ax4.hist(primary[index], bins=50)
    ax4.set_title('primary')

    fig.tight_layout()
    
    plt.savefig(eventfolder+"/trigger_stats.png")
    print("PNG saved:" + eventfolder+"/trigger_stats.png")
    logger.info("PNG saved:" + eventfolder+"/trigger_stats.png")
   
    plt.show()    



    #====== end of run =======      
    logger.info("Done within "+str(time.clock()) +"s")    
    
    
    
