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

from radio_simus.computevoltage import get_voltage, compute_antennaresponse
from radio_simus.signal_processing import run
from radio_simus.in_out import _table_voltage,_load_to_array
#from radio_simus.in_out import inputfromtxt, _get_positions_coreas, inputfromtxt_coreas, load_trace_to_table


if __name__ == '__main__':

    if ( len(sys.argv)<1 ):
        print("""
        Example how to read in an electric field apply antenna response ("voltages") and/or the full chain ("voltages_full")
        
        Usage for full antenna array:
        usage: python example_voltage_hdf5.py <path to event folder> 
        example: python eexample_voltage_hdf5.py ./ coreas
        
        """)
        sys.exit(0)



    ###########################

    directory = sys.argv[1]

    for file in glob.glob(directory+"*.hdf5"):
        print("\n")
        print(file)
        
        ##load info from hdf5 file
        path_hdf5=file
        efield1, time_unit, efield_unit, shower, position, slopes = _load_to_array(path_hdf5, content="efield")
        ## or convert from existing table to numpy array
        #efield1=np.array([a['Time'], a['Ex'], a['Ey'], a['Ez']]).T
        efield=efield1.T
                
        ## apply only antenna response
        #voltage = compute_antennaresponse(efield1, shower['zenith'], shower['azimuth'], alpha=slopes[0], beta=slopes[1] )
        #shower.update({'voltage': 'antennaresponse'})
        #keyword='voltages'
                
        ## apply full chain  --- ToDo make antennaresponse optional
        voltage = run(efield, shower['zenith'], shower['azimuth'], slopes[0], slopes[1], False)
        shower.update({'voltage': ('antennaresponse', 'noise', 'filter', 'digitise')})
        keyword='voltages_full'
                            
        # load voltage array to table and store in same hdf5 file
        volt_table = _table_voltage(voltage, pos=position, slopes=slopes ,info=shower )
        volt_table.write(path_hdf5, path=keyword, format="hdf5", append=True, compression=True,serialize_meta=True) #
