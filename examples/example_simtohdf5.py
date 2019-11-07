################################
#### by A. Zilles:

### As soon as we have decided on 
#           * whether we want to have single file for eachh antenna or one for each event
#           * whether we want to use astropy table or data groups
#   this script can be extremely cleaned-up.

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
#import radio_simus 
#This also loads the submodule in_out, and makes it available without its package prefix
from radio_simus.in_out import inputfromtxt, _get_positions_coreas, inputfromtxt_coreas, load_trace_to_table
from radio_simus.__init__ import Vrms


if __name__ == '__main__':

    if ( len(sys.argv)<2 ):
        print("""
        Script loops over event folders and produces hdf5 files per antenna, containing efield and voltage trace
        
        Example how to read in an electric field or voltage trace, load it to a table, write it into a hdf5 file (and read the file.)
        -- loops over complete event sets
        -- define ALL in script to write all antennas in one hdf5 file and event info in an own column
        -- define SINGLE in script to write single antennas in one hdf5 file each and event info as meta for each antenna
        -- produces hdf5 files for zhaires and coreas simulations
        -- (can compute the full chain for the voltage traces and save them im hdf5 file)
        -- currently only applies antenna response, calculates trigger and p2p, and stores them in the hadf5 file 
        -- example how to read in tables from hdf5 files
        
        Usage for full antenna array:
        usage: python example_simtohdf5.py <path to event folders> <zhaires/coreas>
        example: python example_simtohdf5.py ./ coreas
        
        ATTENTION: 
        -- adapt the paths given in the config-file so that eg Vrms (for threshold) can be read-in correctly
        -- To be specified SINGLE/ALL True or False depending on favorite saving mode, same size file in total
        -- Hardcoded threshold value in muV, read-in from config file
        
        Note: 
        --- ToDo: adopt example to also add voltage traces from txt files
        --- ToDo: save units of parameters, astropy units or unyt
        --- ToDo: implement calculation of alpha and beta, so far read in from file or set to (0,0) [in_out.py]
        --- ToDo: fix the problem of array size in computevoltage which leads to not-produced traces at the moment
        """)
        sys.exit(0)
    
    
    
    
    # path to folder containing traces
    eventfolder = sys.argv[1]
    # coreas or zhaires -- treatment differently
    simus = str(sys.argv[2])
   
    SINGLE=True # save each antenna as single antenna file
    ALL=False # all antenna in one file
   
    # for triggering
    threshold_aggr = 3*Vrms #muV
    threshold_cons = 5*Vrms #muV  
    
    # set option on what done in this script, more options in the script
    # here, we switch on teh computation of the voltage response at teh storage of the voltage trace in the hdf5 file
    # + trigger info + p2p value (at least for now)  -- not done (digitise, filter etc, can be done as well)
    voltage_compute=True
    if voltage_compute:
        from radio_simus.computevoltage import get_voltage, compute_antennaresponse
        from radio_simus.signal_processing import run
        from radio_simus.in_out import _table_voltage
        
        
        
    #----------------------------------------------------------------------   
    # General info statement for later identification
    logger.info("general: Eventset =" +str(eventfolder)+ ", Simulation = "  +str(simus)+ ", SINGLE/ALL = "  +str(SINGLE)+ "/"  +str(ALL)+  ", Threshold = "  +str(threshold_aggr)+ "(aggr), "  +str(threshold_cons)+ "(cons)")
    print("\nGENERAL: Eventset =" +str(eventfolder)+ ", Simulation = "  +str(simus)+ ", SINGLE/ALL = "  +str(SINGLE)+ "/"  +str(ALL)+  ", Threshold = "  +str(threshold_aggr)+ " (aggr), "  +str(threshold_cons)+ " (cons)\n")
    #----------------------------------------------------------------------   
   
   
    # loop over eventfolder
    for path in glob.glob(eventfolder+"/*/"):
        
      # check whether it is a folder
      if os.path.isdir(path): 
   
        showerID=str(path).split("/")[-2] # should be equivalent to folder name
        print(" --- Event : ", showerID, ", in ", path)
        logger.debug(" --- Event : "+ str(showerID)+ ", in "+ str(path) )
        
        task=None
        core=None    
        
        
        ##########################
        ### LOADING EVENT INFO ###
        ##########################

        
        if simus == 'zhaires':
            ####################################### NOTE zhaires
            # Get the antenna positions from file
            positions = np.loadtxt(path+"antpos.dat")
            ID_ant = []
            slopes = []
            # TODO adopt reading in positions, ID_ant and slopes to coreas style - read in from SIM*info    
            #posfile = path +'SIM'+str(showerID)+'.info'
            #positions, ID_ant, slopes = _get_positions_coreas(posfile)
            ##print(positions, ID_ant, slopes)
                
            # Get shower info
            inputfile = path+showerID+'.inp'
            #inputfile = path+"/inp/"+showerID+'.inp'
            #print("Check inputfile path: ", inputfile)
            try:
                zen,azim,energy,injh,primarytype,core,task = inputfromtxt(inputfile)
            except:
                print("no TASK, no CORE")
                inputfile = path+showerID+'.inp'
                zen,azim,energy,injh,primarytype = inputfromtxt(inputfile)
            
            # correction of shower core
            positions = positions + np.array([core[0], core[1], 0.])

            ending_e = "a*.trace"
                

        if  simus == 'coreas':
            #posfile = path +'SIM'+str(showerID)+'.list' # contains not alpha and beta
            posfile = path +'SIM'+str(showerID)+'.info' # contains original ant ID , positions , alpha and beta
            positions, ID_ant, slopes = _get_positions_coreas(posfile)
            #print(positions, ID_ant, slopes)
            
            inputfile = path+'/inp/SIM'+showerID+'.inp'
            zen,azim,energy,injh,primarytype,core,task = inputfromtxt_coreas(inputfile)
            
            # correction of shower core
            positions = positions + np.array([core[0]/u.m, core[1]/u.m, 0.*u.m/u.m])

            #redefinition of path to traces
            path=path+'/SIM'+showerID+'_coreas/'
            
            ending_e = "raw_a*.dat"
            
    
        #----------------------------------------------------------------------   

    
        ########################
        # load shower info from inp file via dictionary
        ########################
        shower = {
                "ID" : showerID,               # shower ID, number of simulation
                "primary" : primarytype,        # primary (electron, pion)
                "energy" : energy,               # eV
                "zenith" : zen,               # deg (GRAND frame)
                "azimuth" : azim,                # deg (GRAND frame)
                "injection_height" : injh,    # m (injection height in the local coordinate system)
                "task" : task,    # Identification
                "core" : core,    # m, numpy array, core position
                "simulation" : simus # coreas or zhaires
                }
        ####################################
        print("shower", shower)
        logger.info("Shower summary: " + str(shower))
        
        
        #shower.write(name_all, path='event', format="hdf5", append=True,  compression=True,serialize_meta=True) 
        #positions.write(name_all, path='positions', format="hdf5", append=True,  compression=True,serialize_meta=True) 
        #slopes.write(name_all, path='slopes', format="hdf5", append=True,  compression=True,serialize_meta=True) 
        #ID_ant.write(name_all, path='IDs', format="hdf5", append=True,  compression=True,serialize_meta=True)
        
        if ALL:
            name_all = path+'/event_'+showerID+'.hdf5'
            ##hf = h5py.File(name_all, 'w')
            #hf = h5py.File(name_all, 'w')
            #hf.create_dataset('positions', data=positions)
            #hf.create_dataset('slopes', data=slopes)
            ##hf.create_dataset('ID_ant', data=np.asarray(ID_ant))
            ##hf.create_dataset('shower', data=shower)
            ##dset = hf.create_dataset("shower", shower) 
            #hf.close()
    
            a1 = Column(data=np.array(ID_ant), name='ant_ID')

            b1 = Column(data=positions.T[0], unit=u.m, name='pos_x')
            c1 = Column(data=positions.T[1], unit=u.m, name='pos_y')
            d1 = Column(data=positions.T[2], unit=u.m, name='pos_z')  #u.eV, u.deg
            e1 = Column(data=slopes.T[0], unit=u.deg, name='alpha')
            f1 = Column(data=slopes.T[1], unit=u.deg, name='beta')  #u.eV, u.deg    
            event_info = Table(data=(a1,b1,c1,d1,e1,f1,), meta=shower) 
            event_info.write(name_all, path='event', format="hdf5", append=True,  compression=True, serialize_meta=True)



        for ant in tqdm.tqdm(glob.glob(path+ending_e)):
            #print("\n Read in data from ", ant)
            
            name = None
            if simus == 'zhaires':
                ant_number = int(ant.split('/')[-1].split('.trace')[0].split('a')[-1]) # index in selected antenna list
                ID = ant_number # original ant ID
                # TODO adopt to read in SIM*info file, reprocude ant ID in orginal list
                #ant_number = int(ant.split('/')[-1].split('.trace')[0].split('a')[-1])
                # ID = ID_ant.index(ant_number) # index of antenna in positions
                # define path for storage of hdf5 files
                if SINGLE:
                    name = path+'/table_'+str(ID)+'.hdf5' # to adopt to coreas style for now 
                    #print("-- ant_ID: "+ str(ID) + "-- Table saved as: " +str(name))
                    logger.debug("-- ant_ID: "+ str(ID) + "-- Table saved as: " +str(name))

            if  simus == 'coreas':
                base=os.path.basename(ant)
                #coreas
                #ID == ant_number for coreas
                ID=(os.path.splitext(base)[0]).split("_a")[1] # remove raw , original ant ID
                ant_number = ID_ant.index(ID) # index of antenna in positions list
                # define path for storage of hdf5 files
                if SINGLE:
                    name = path+'/../table_'+str(ID)+'.hdf5'
                    #print("-- ant_ID: "+ str(ID) + "-- Table saved as: " +str(name))
                    logger.debug("-- ant_ID: "+ str(ID) + "-- Table saved as: " +str(name))
                    
            ##### read-in output of simulations        
            if SINGLE:
                # read in trace from file and store as astropy table, saved as hdf5 file (optional), conversion of units done internally
                # a is a astropy table, saved as <name> as hdf5-file
                a= load_trace_to_table(path=ant, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist(),  info=shower, content="e", simus=simus, save=name) 
                
                ## write astropy table to hdf5 file
                ###a.write(name, path=name, overwrite=True, compression=True, serialize_meta=True) #append=True, 
                #a.write(name, path='efield', format="hdf5",  serialize_meta=True)
            
                ####### From here on only additional feature and nice-to-know
                ## play around
                #print(a.info)
                #print(a['Ex'])
                #print(a.meta)
                #print(type(a))
                #print(a)
            
            if ALL:
                # create a group per antenna  --- save ID, position, slope as attributes of group 
                #g1 = hf.create_group(str(ant_number))
                
                a= load_trace_to_table(path=ant, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist(),  info={}, content="e", simus=simus, save=name_all, ant="/"+str(ID)+"/")
                
                #g1.create_dataset('efield', data = a,compression="gzip", compression_opts=9)
                
                #efield = np.loadtxt(ant)
                #if simus=="coreas": 
                    #efield.T[0]*=1e9 # s to ns
                    ### coreas cgs to SI, V/m to muV/m
                    #efield.T[1]*=2.99792458e4* 1.e6 
                    #efield.T[2]*=2.99792458e4* 1.e6 
                    #efield.T[3]*=2.99792458e4* 1.e6
                #g1.create_dataset('efield', data = efield,compression="gzip", compression_opts=9)
            

            ####################################################################

            ############
            ## PLAYGROUND
            ############
            
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
                
                
            ########### example how to add voltage from text file as a column    
            voltage=False
            if voltage:
                ##### Hopefully not needed any more if voltage traces are not stored as txt files in future
                print("Adding voltages")
                
                ### ATTENTION currently electric field added - adjust <ant=path+ending_e>
                print("WARNING: adopt path to voltage trace")
                logger.warning("Read in voltage trace from file : adopt path to voltage trace")
                
                if SINGLE:
                    # read in trace from file and store as astropy table - can be substituted by computevoltage operation
                    b= load_trace_to_table(path=ant, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist(), info=shower, content="v") # read from text files
        
                    ## stack electric field and voltage traces
                    #from astropy.table import hstack, then write c to a file
                    #c = hstack([a, b])
                
                    # Write to tables in hdf5 file of the efield
                    b.write(name, path='voltages', append=True, serialize_meta=True) #append=True -- NOTE: Do I need that
                
                
                if ALL:
                    # read in trace from file and store as astropy table - can be substituted by computevoltage operation
                    b= load_trace_to_table(path=ant, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist(), info=None, content="v",simus=simus, save=name_all, ant="/"+str(ID)+"/") # read from text files
                    #g1.create_dataset('voltages', data = b,compression="gzip", compression_opts=9)
            
            ########### example VOLTAGE COMPUTATION and add to same hdf5 file
            #voltage_compute=True
            if voltage_compute:
                
                ##load info from hdf5 file
                #path_hdf5=name
                #efield1, time_unit, efield_unit, shower, position = _load_to_array(path_hdf5, content="efield")
                # or convert from existing table to numpy array
                efield1=np.array([a['Time'], a['Ex'], a['Ey'], a['Ez']]).T
                
                
                try:
                    ## apply only antenna response
                    voltage = compute_antennaresponse(efield1, shower['zenith'], shower['azimuth'], alpha=slopes[ant_number,0], beta=slopes[ant_number,1] )

                    ## NOTE apply full chain, not only antenna resonse: add noise, filter, digitise
                    #voltage = run(efield1, shower['zenith'], shower['azimuth'], 0, 0, False) # alpha = 0, beta = 0
                    
                    ### add some info on P2P and  TRIGGER if wanted: trigger on any component, or x-y combined
                    #NOTE: need to cross-check that trace is also in muV, but should be by definition
                    from radio_simus.signal_treatment import p2p, _trigger
                    # peak-to-peak values: x, y, z, xy-, all-combined
                    p2p_values = p2p(voltage)
                    # trigger info: trigger = [thr_aggr, any_aggr, xz_aggr, thr_cons, any_cons, xy_cons]
                    trigger =  [threshold_aggr, _trigger(p2p_values, 'any', threshold_aggr/(u.u*u.V)), _trigger(p2p_values, 'xy', threshold_aggr/(u.u*u.V)), threshold_cons, _trigger(p2p_values, 'any', threshold_cons/(u.u*u.V)), _trigger(p2p_values, 'xy', threshold_cons/(u.u*u.V))]
                    
                    # Update info
                    shower.update({'voltage': 'antennaresponse', 'trigger': trigger, 'p2p': list(p2p_values) })
                    #shower.update({'voltage': ('antennaresponse', 'noise', 'filter', 'digitise')})

                    if SINGLE:   
                        # load voltage array to table and store in same hdf5 file
                        volt_table = _table_voltage(voltage, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist() ,info=shower )
                        volt_table.write(name, path='voltages', format="hdf5", append=True, serialize_meta=True) #compression=True,

                            
                    if ALL:
                        #voltage = compute_antennaresponse(efield, shower['zenith'], shower['azimuth'], alpha=slopes[ant_number,0], beta=slopes[ant_number,1] )
                        #g1.create_dataset('voltages', data = volt_table, compression="gzip", compression_opts=9)
                        
                        #v_info = {'voltage': 'antennaresponse'}
                        #v_info = {'voltage': ('antennaresponse', 'noise', 'filter', 'digitise')}
                        volt_table = _table_voltage(voltage, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist() ,info=shower) #v_info )

                        volt_table.write(name_all, path="/"+str(ID)+"/"+'voltages', format="hdf5", append=True, compression=True,serialize_meta=True) 
                        
                except: 
                        print("====== ATTENTION: ValueError raised for a"+str(ant_number) + " --- check computevoltage =======")
                        logger.error("ValueError raised for a"+str(ant_number) + "--- check computevoltage")



            ######## just testing part and examples how to use astropy tables
            EXAMPLE=False
            if EXAMPLE:
                if SINGLE:
                    # read in hdf5 file 
                    f=Table.read(name, path="efield")
                    #g=Table.read(name, path="voltages")
                    print(f)
                    #print(g)

                    #b=f['Ex']*f['Ex']
                    #print(b)
                    print(f.meta, f.info)
                    #print(f.meta['zenith'], f['Time'].unit)
                    #print(g)
                    
                    ### Just examples how output could handled 
                    #summe=f['Ex']+f['Ey']
                    #print(summe[-2])
                if ALL:
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
