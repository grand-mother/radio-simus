################################
#### by A. Zilles:

### this scripts reads in a raw simulations and saves it as a hdf5 file
### its also include a way to compute the voltage trace
### traces and infos saved in one singles hdf5 file, as astropy tables

################################
#!/usr/bin/env python


import os
from os.path import split, join, realpath
import numpy as np
import sys
import glob

import logging
logging.basicConfig(filename="example_simtohdf5.log", level=logging.INFO)
logger = logging.getLogger('Main')

import tqdm

import astropy
from astropy.table import Table, Column
from astropy.table import hstack
from astropy import units as u  
import h5py


# Expand the PYTHONPATH
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))

# load config file first
import radio_simus 
radio_simus.load_config('./test.config')
Vrms = radio_simus.config.Vrms

from radio_simus.io_utils import load_trace_to_table, load_eventinfo_tohdf


if __name__ == '__main__':

    if ( len(sys.argv)<2 ):
        print("""
        Script loops over event folders and produces hdf5 files per antenna, containing efield and voltage trace
        Note that this is part of the restructering phase for the software, large fraction of this scripts can be substituted by module which hide the action for the end-user. But first, we have to agree one a common approach...
        
        Example how to read in an electric field or voltage trace, load it to a table, write it into a hdf5 file (and read the file.)
        -- loops over complete event sets
        -- DEFAULT: define ALL in script to write all antennas in one hdf5 file and event info in an own column
        ##-- define SINGLE in script to write single antennas in one hdf5 file each and event info as meta for each antenna
        -- produces hdf5 files for zhaires and coreas simulations
        -- (can compute the full chain for the voltage traces and save them im hdf5 file)
        -- currently only applies antenna response, calculates trigger and p2p, and stores them in the hadf5 file 
        -- example how to read in tables from hdf5 files
        
        Usage for full antenna array:
        usage: python example_simtohdf5.py <path to event folders> <zhaires/coreas>
        example: python example_simtohdf5.py ./ coreas
        
        ATTENTION: 
        -- adapt the paths given in the config-file so that eg Vrms (for threshold) can be read-in correctly
        -- only all antennas in one hdf5 file possible
        -- Hardcoded ; in voltage_compute option: list of appllied processing modules: processing_info={..}
        
        Note: 
        --- ToDo: save units of parameters, astropy units or unyt
        --- ToDo: implement calculation of alpha and beta, so far read in from file or set to (0,0) [io_utils.py]
        --- ToDo: fix the problem of array size in computevoltage which leads to not-produced traces at the moment
        """)
        sys.exit(0)
    
    
    
    
    # path to folder containing traces
    eventfolder = sys.argv[1]
    # coreas or zhaires -- treatment differently
    simus = str(sys.argv[2])
   
   
    # for triggering
    threshold_aggr = 3*Vrms #muV
    threshold_cons = 5*Vrms #muV  
    
    # set option on what done in this script, more options in the script
    # here, we switch on teh computation of the voltage response at teh storage of the voltage trace in the hdf5 file
    # + trigger info + p2p value (at least for now)  -- not done (digitise, filter etc, can be done as well)
    voltage_compute=True
    if voltage_compute:
        from radio_simus.computevoltage import get_voltage, compute_antennaresponse
        from radio_simus.signal_processing import standard_processing
        from radio_simus.io_utils import _table_voltage
        
        
        
    #----------------------------------------------------------------------   
    # General info statement for later identification
    logger.info("general: Eventset =" +str(eventfolder)+ ", Simulation = "  +str(simus)#+ ", SINGLE/ALL = "  +str(SINGLE)+ "/"  +str(ALL)
                +", Threshold = "  +str(threshold_aggr)+ "(aggr), "  +str(threshold_cons)+ "(cons)")
    print("\nGENERAL: Eventset =" +str(eventfolder)+ ", Simulation = "  +str(simus) #+ ", SINGLE/ALL = "  +str(SINGLE)+ "/"  +str(ALL)
          +  ", Threshold = "  +str(threshold_aggr)+ " (aggr), "  +str(threshold_cons)+ " (cons)\n")
    #----------------------------------------------------------------------   
   
   
    # loop over eventfolder
    for path in glob.glob(eventfolder+"/*/"):
        
      # check whether it is a folder
      if os.path.isdir(path): 
   
        showerID=str(path).split("/")[-2] # should be equivalent to folder name
        print(" --- Event : ", showerID, ", in ", path)
        logger.debug(" --- Event : "+ str(showerID)+ ", in "+ str(path) )
           
        ### Name of hdf5 file for storage
        name_all = path+'/../event_'+showerID+'.hdf5'
        

        
        ##########################
        ### LOADING EVENT INFO ###
        ##########################
        
        shower, ID_ant, positions, slopes = load_eventinfo_tohdf(path, showerID, simus, name_all)

        
        if simus == "zhaires":
            #redefinition of path to traces needed.. depends on folder structure
            #
            # trace ending
            ending_e = "a*.trace"

        if  simus == 'coreas':
            #redefinition of path to traces
            path=path+'/SIM'+showerID+'_coreas/'
            # trace ending
            ending_e = "raw_a*.dat"
            
        ########################
        ### EXAMPLE ANALYSIS ###
        ########################
        p2p_values=[]
        trigger=[]    
            
        # loop over existing single antenna files as raw output from simulations
        for ant in tqdm.tqdm(glob.glob(path+ending_e)):
            #print("\n Read in data from ", ant)

            #### Get correct antenna number==index in numpy arrays from underlying ID of antenna
            ## ID = antenna identification, ant_number == index in array
            if simus == 'zhaires':
                ant_number = int(ant.split('/')[-1].split('.trace')[0].split('a')[-1]) # index in selected antenna list
                ID = ant_number # original ant ID
                # TODO adopt to read in SIM*info file, reprocude ant ID in orginal list
                #ant_number = int(ant.split('/')[-1].split('.trace')[0].split('a')[-1])
                # ID = ID_ant.index(ant_number) # index of antenna in positions
                # define path for storage of hdf5 files
                

            if  simus == 'coreas':
                base=os.path.basename(ant)
                #coreas
                #ID == ant_number for coreas
                ID=(os.path.splitext(base)[0]).split("_a")[1] # remove raw , original ant ID
                ant_number = ID_ant.index(ID) # index of antenna in positions list
                # define path for storage of hdf5 files
                #print(ant_number, ID)

                    
            ##### read-in output of simulations and save as hdf5 file     
            a= load_trace_to_table(path_raw=ant, pos=positions[ant_number].value.tolist(), slopes=slopes[ant_number].value.tolist(), info={}, content="e", simus=simus, save=name_all, ant="/"+str(ID)+"/")

            

            ###################################################################
            ########### example VOLTAGE COMPUTATION and add to same hdf5 file

            if voltage_compute:
                
                ##load info from hdf5 file
                #efield1, time_unit, efield_unit, shower, position = _load_to_array(path_hdf5, content="efield", ant="/"+str(ID)+"/")
                
                # or convert from existing table to numpy array
                efield1=np.array([a['Time'], a['Ex'], a['Ey'], a['Ez']]).T

                try:
                    ## apply only antenna response -- the standard applying the antenna response
                    processing_info={'voltage': 'antennaresponse'} # HARDCODED
                    voltage = compute_antennaresponse(efield1, shower['zenith'].value, shower['azimuth'].value, 
                                                      alpha=slopes[ant_number,0].value, beta=slopes[ant_number,1].value)
                    

                    ## NOTE apply full chain, not only antenna resonse: add noise, filter, digitise
                    #processing_info={'voltage': ('antennaresponse', 'noise', 'filter', 'digitise')}
                    #voltage=standard_processing(efield1, shower['zenith'].value, shower['azimuth'].value, 
                                                #alpha=slopes[ant_number,0].value, beta=slopes[ant_number,1].value,
                                                #processing=processing_info["voltage"], DISPLAY=0)
                    
                    # convert to astropy table and save in hdf5 files
                    volt_table = _table_voltage(voltage, pos=positions[ant_number].value.tolist(), slopes=slopes[ant_number].value.tolist()     
                                                ,info=processing_info, save=name_all, ant="/"+str(ID)+"/")
                        
                except: 
                    print("====== ATTENTION: ValueError raised for a"+str(ant_number) + " --- check computevoltage =======")
                    logger.error("ValueError raised for a"+str(ant_number) + "--- check computevoltage")
                
                
                ##############################
                #### DO A FIRST ANALYSIS #####
                ##############################
                
                try:
                    ### add some info on P2P and  TRIGGER if wanted: trigger on any component, or x-y combined
                    #NOTE: need to cross-check that trace is also in muV, but should be by definition
                    from radio_simus.signal_treatment import p2p, _trigger
                    
                    # peak-to-peak values: x, y, z, xy-, all-combined
                    p2p_values.append(p2p(voltage))

                    # trigger info: trigger = [thr_aggr, any_aggr, xz_aggr, thr_cons, any_cons, xy_cons]
                    trigger.append( (_trigger(p2p(voltage), 'any', threshold_aggr/(u.u*u.V)), _trigger(p2p(voltage), 'xy', threshold_aggr/(u.u*u.V)),  _trigger(p2p(voltage), 'any', threshold_cons/(u.u*u.V)), _trigger(p2p(voltage), 'xy', threshold_cons/(u.u*u.V))) )
        
                except: 
                    print("====== ATTENTION: ValueError raised for a"+str(ant_number) + " --- no first analysis performed =======")
                    logger.error("ValueError raised for a"+str(ant_number) + "--- check computevoltage")
        
        
        ## Quit loop 
        p2p_values = np.array(p2p_values)
        trigger = np.array(trigger)
            
        a2 = Column(data=np.array(ID_ant), name='ant_ID')
                        
        b2 = Column(data=p2p_values.T[0], unit=u.u*u.V, name='p2p_x')  
        c2 = Column(data=p2p_values.T[1], unit=u.u*u.V, name='p2p_y') 
        d2 = Column(data=p2p_values.T[2], unit=u.u*u.V, name='p2p_z') 
        e2 = Column(data=p2p_values.T[3], unit=u.u*u.V, name='p2p_xy') 
                        
        f2 = Column(data=trigger.T[0],  name='trigger_aggr_any')
        g2 = Column(data=trigger.T[1],  name='trigger_aggr_xy')
        h2 = Column(data=trigger.T[2],  name='trigger_cons_any')  
        i2 = Column(data=trigger.T[3],  name='trigger_cons_xy')
                        
        thres_info = {"threshold_aggr": threshold_aggr, "threshold_cons": threshold_cons }
        analysis_info = Table(data=(a2,b2,c2,d2,e2,f2,g2,h2,i2,), meta=thres_info) 
        analysis_info.write(name_all, path='analysis', format="hdf5", append=True,  compression=True, serialize_meta=True)


#----------------------------------------------------------------------   



''' must be in the loop
            ############
            ## PLAYGROUND
            ############
                
                
            ############ example how to add voltage from text file as a column    
            #voltage=False
            #if voltage:
                ###### Hopefully not needed any more if voltage traces are not stored as txt files in future
                #print("Adding voltages")
                
                #### ATTENTION currently electric field added - adjust <ant=path+ending_e>
                #print("WARNING: adopt path to voltage trace")
                #logger.warning("Read in voltage trace from file : adopt path to voltage trace")

                ## read in trace from file and store as astropy table - can be substituted by computevoltage operation
                ##b= load_trace_to_table(path=ant, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist(), info=None, content="v",simus=simus, save=name_all, ant="/"+str(ID)+"/") # read from text files
                #b= load_trace_to_table(path=ant, pos=positions[ant_number].value.tolist(), slopes=slopes[ant_number].value.tolist(), info=None, content="v",simus=simus, save=name_all, ant="/"+str(ID)+"/") # read from text files
                ##g1.create_dataset('voltages', data = b,compression="gzip", compression_opts=9)



            ########## Plotting a astropy table
            DISPLAY=False
            if DISPLAY:
                import matplotlib.pyplot as plt

                plt.figure(1,  facecolor='w', edgecolor='k')
                plt.plot(a['Time'],a['Ex'], label="Ex")
                plt.plot(a['Time'],a['Ey'], label="Ey")
                plt.plot(a['Time'],a['Ez'], label="Ez")

                plt.xlabel('Time ('+str(a['Time'].unit)+')')
                plt.ylabel('Electric field ('+str(a['Ex'].unit)+')')
                plt.legend(loc='best')
                
                plt.show()
                #plt.savefig('test.png', bbox_inches='tight')

            ######## just testing part and examples how to use astropy tables
            EXAMPLE=False
            if EXAMPLE:
                #with h5py.File(name_all, 'r') as f:
                    ##group = f[str(ant_number)]
                    ##dataset = group['efield']
                    ##print(dataset)
                    ##print('List of items in the base directory:', f.items(), f.attrs['ID'])
                    #info = f.get("info")
                    #gp1 = f.get(str(ant_number))
                    #gp2 = f.get(str(ant_number)+'/efield')
                    ##print(u'\nlist of items in the group 1:', gp1.items())
                    ##print(u'\nlist of items in the subgroup:', gp2.items())
                    #arr1 = np.array(gp1.get('efield'))
                    #print(info, arr1)
                        
                f=Table.read(name_all, path="/"+str(ID)+"/efield")
                g=Table.read(name_all, path="/event")
                print(f)
                print(f.meta, f.info, g.meta, g.info)
'''
