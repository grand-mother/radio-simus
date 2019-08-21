################################
#### by A. Zilles
################################
#!/usr/bin/env python


import os
from os.path import split, join, realpath
import numpy as np
import sys
import glob

#from in_out import _get_positions_coreas, inputfromtxt_coreas, load_trace_to_table
# Expand the PYTHONPATH and import the radiomorphing package #NOTE: this would be on the shared disc
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
#import radio_simus 
#This also loads the submodule in_out, and makes it available without its package prefix
from radio_simus.in_out import inputfromtxt, _get_positions_coreas, inputfromtxt_coreas, load_trace_to_table


if __name__ == '__main__':

    if ( len(sys.argv)<1 ):
        print("""
        Example how to read in an electric field or voltage trace, load it to a table, write it into a hdf5 file (and read the file.)
        -- produces hdf5 files for zhaires and coreas simulations
        -- can compute the full chain for the voltage traces and save them im hdf5 file
        -- can only only apply antenna response
        -- example hoe to read in tables from hdf5 files
        
        Usage for full antenna array:
        usage: python example_simtohdf5.py <path to event folder> <zhaires/coreas>
        example: python example_simtohdf5.py ./ coreas
        
        ATTENTION: To be specified SINGLE/ALL True or False depending on favorite saving mode
        
        Note: adopt example to also add voltage traces from txt files
        """)
        sys.exit(0)
    
    # path to folder containing traces
    path = sys.argv[1]
    # coreas or zhaires -- treatment differently
    simus = str(sys.argv[2])
   
    SINGLE=None # save each antenna as single antenna file
    ALL=True # all antenna in one file
   
   
   #----------------------------------------------------------------------   
   
   
    showerID=str(path).split("/")[-2] # should be equivalent to folder name
    print(" --- Event : ", showerID, ", in ", path)
    
    task=None
    core=None    
    if simus == 'zhaires':
        ####################################### NOTE zhaires
        # Get the antenna positions from file
        positions = np.loadtxt(path+"antpos.dat")
            
        # Get shower info
        inputfile = path+showerID+'.inp'
        #inputfile = path+"/inp/"+showerID+'.inp'
        print("Check inputfile path: ", inputfile)
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
        posfile = path +'SIM'+str(showerID)+'.info'
        positions, ID_ant, slopes = _get_positions_coreas(posfile)
        #print(positions, ID_ant, slopes)
        
        inputfile = path+'/inp/SIM'+showerID+'.inp'
        zen,azim,energy,injh,primarytype,core,task = inputfromtxt_coreas(inputfile)
        
        # correction of shower core
        positions = positions + np.array([core[0], core[1], 0.])

        #redefinition of path to traces
        path=path+'/SIM'+showerID+'_coreas/'
        
        ending_e = "raw_a*.dat"
        
 
 
 
    from astropy import units as u
    print("ToDo: save units of parameters, astropy units or unyt")
    print("ToDo: save shower parameters on highest level, together with availble antenn positions")
    
    ########################
    # load shower info from inp file via dictionary
    ########################
    shower = {
            "ID" : showerID,               # shower ID, number of simulation
            "primary" : primarytype,        # primary (electron, pion)
            "energy" : energy,               # EeV
            "zenith" : zen,               # deg (GRAND frame)
            "azimuth" : azim,                # deg (GRAND frame)
            "injection_height" : injh,    # m (injection height in the local coordinate system)
            "task" : task,    # Identification
            "core" : core.tolist(),    # m, numpy array, core position
            "simulation" : simus # coreas or zhaires
            }
    ####################################
    print("shower", shower)
    
    
    shower.write(name_all, path='event', format="hdf5", append=True,  compression=True,serialize_meta=True) 
    
    import astropy
    from astropy.table import Table
    from astropy.table import hstack
    import h5py

    
    if ALL:
        name_all = path+'/event_'+showerID+'.hdf5'
        #hf = h5py.File(name_all, 'w')
        
   

    for ant in glob.glob(path+ending_e):
        print("\n Read in data from ", ant)
        
        name = None
        if simus == 'zhaires':
            ant_number = int(ant.split('/')[-1].split('.trace')[0].split('a')[-1])
            # define path for storage of hdf5 files
            if SINGLE:
                name = path+'/table_a'+str(ant_number)+'.hdf5' # to adopt to coreas style for now 
                print("Table saved as: ", name)

        if  simus == 'coreas':
            base=os.path.basename(ant)
            # coreas
            ID=(os.path.splitext(base)[0]).split("_a")[1] # remove raw 
            ant_number = ID_ant.index(ID)
            # define path for storage of hdf5 files
            if SINGLE:
                name = path+'/../table_'+str(ID)+'.hdf5'
                print("Table saved as: ", name)
                
         ##### read-in output of simulations        
        if SINGLE:
            # read in trace from file and store as astropy table, saved as hdf5 file (optional)
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
            
            a= load_trace_to_table(path=ant, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist(),  info=shower, content="e", simus=simus, save=name_all, ant="/"+str(ant_number)+"/")
            
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
        
        ########## Plotting a table
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
                b= load_trace_to_table(path=ant, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist(), info=shower, content="v",simus=simus, save=name_all, ant="/"+str(ant_number)+"/") # read from text files
                #g1.create_dataset('voltages', data = b,compression="gzip", compression_opts=9)
        
        ########### example VOLTAGE COMPUTATION and add to same hdf5 file
        voltage_compute=True
        if voltage_compute:
            from radio_simus.computevoltage import get_voltage, compute_antennaresponse
            from radio_simus.signal_processing import run
            from radio_simus.in_out import _table_voltage
            
            
            ##load info from hdf5 file
            #path_hdf5=name
            #efield1, time_unit, efield_unit, shower, position = _load_to_array(path_hdf5, content="efield")
            # or convert from existing table to numpy array
            efield1=np.array([a['Time'], a['Ex'], a['Ey'], a['Ez']]).T
                
                
            ## apply only antenna response
            voltage = compute_antennaresponse(efield1, shower['zenith'], shower['azimuth'], alpha=slopes[ant_number,0], beta=slopes[ant_number,1] )
            shower.update({'voltage': 'antennaresponse'})
                
            ## apply full chain
            #voltage = run(efield1, shower['zenith'], shower['azimuth'], 0, 0, False)
            #shower.update({'voltage': ('antennaresponse', 'noise', 'filter', 'digitise')})
                
            # load voltage array to table and store in same hdf5 file
            volt_table = _table_voltage(voltage, pos=positions[ant_number].tolist(), slopes=slopes[ant_number].tolist() ,info=shower )
            
            
            if SINGLE:             
                volt_table.write(name, path='voltages', format="hdf5", append=True, serialize_meta=True) #compression=True,

            if ALL:
                #voltage = compute_antennaresponse(efield, shower['zenith'], shower['azimuth'], alpha=slopes[ant_number,0], beta=slopes[ant_number,1] )
                #g1.create_dataset('voltages', data = volt_table, compression="gzip", compression_opts=9)
                
                volt_table.write(name_all, path="/"+str(ant_number)+"/"+'voltages', format="hdf5", append=True, compression=True,serialize_meta=True) 




        ######## just testing part and examples how to use astropy tables
        EXAMPLE=True
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
                    
                f=Table.read(name_all, path="/"+str(ant_number)+"/efield")
                print(f)
                print(f.meta, f.info)
