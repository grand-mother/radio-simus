import sys
from sys import argv
import os
import glob

import numpy as np
from numpy import *


#from astropy import units as u
#muV_m = u.u * u.V / u.m

from astropy.table import Table
from astropy.table import hstack
import h5py

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))

from radio_simus.computevoltage import  compute_antennaresponse
from radio_simus.signal_processing import standard_processing
from radio_simus.io_utils import _table_voltage,_load_to_array, _load_eventinfo_fromhdf, _load_path


if __name__ == '__main__':

    if ( len(sys.argv)<1 ):
        print("""
        Example how to read in an electric field from hdf5 
        apply antenna response ("voltages":'antennaresponse') or the full chain ("voltages": (...) )
        and store the voltage trace in hdf5 file
        
        Usage for full antenna array:
        usage: python example_voltage_hdf5.py <path to event folder> 
        example: python example_voltage_hdf5.py ./ 
        
        Note: signal processing chain hardcoded, see 
            processing_info={'voltage': ('antennaresponse', 'noise', 'filter', 'digitise')}
        
        """)
        sys.exit(0)



    ###########################

    #path to folder with hdf5 event files
    directory = sys.argv[1]
    
    ## choose the steps in processing of electric field
    #processing_info={'voltage': 'antennaresponse'} # calls only compute_antennaresponse
    processing_info={'voltage': ('antennaresponse', 'noise', 'filter', 'digitise')}


    for file in glob.glob(directory+"*.hdf5"): # loop over files in folder
        print("\n")
        print(file)

        ## get information on event from hdf5 file -- in principle not needed here
        shower, ant_ID, positions, slopes = _load_eventinfo_fromhdf(file)
        #print( len(positions[0]) , "antenna positions simulated in Event ", shower["ID"])
        
        ## get information on analysis already performed
        #analysis = _load_path(file, path="/analysis")
        
        f = h5py.File(file, 'a') # open event files
        for ID in f.keys(): #loop over antennas in event
            try:
                print(ID)
                # obtain electric field trace from hdf5 file
                efield, time_unit, efield_unit, ant_position, ant_slopes = _load_to_array(file, content="efield", ant=ID)
                
                # apply voltage treatment
                voltage=standard_processing(efield, shower['zenith'], shower['azimuth'], alpha_sim=ant_slopes[0], beta_sim=ant_slopes[1],
                                    processing=processing_info["voltage"], DISPLAY=0)
                # save trace in hdf5 file
                volt_table = _table_voltage(voltage, pos=ant_position, slopes=ant_slopes ,info=processing_info, 
                                            save=file, ant="/"+str(ID)+"/") 
                
                
            except: # skips all the keys which are not antenna files
                continue
            
                

