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
import radio_simus 
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
        usage: python in_out.py <path to event folder> <zhaires/coreas>
        example: python in_out.py ./ coreas
        
        Note: adopt example to also add voltage traces from txt files
        """)
        sys.exit(0)
    
    # path to folder containing traces
    path = sys.argv[1]
    
    # coreas or zhaires -- treatment differently
    simus = str(sys.argv[2])
   
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
        
        # coorection of shower core
        positions = positions + np.array([core[0], core[1], 0.])

        ending_e = "a*.trace"
               

    if  simus == 'coreas':
        posfile = path +'SIM'+str(showerID)+'.list'
        positions, ID_ant = _get_positions_coreas(posfile)
        
        inputfile = path+'/inp/SIM'+showerID+'.inp'
        zen,azim,energy,injh,primarytype,core,task = inputfromtxt_coreas(inputfile)
        
        # coorection of shower core
        positions = positions + np.array([core[0], core[1], 0.])

        #redefinition of path to traces
        path=path+'/SIM'+showerID+'_coreas/'
        
        ending_e = "raw_a*.dat"
        
        
    ########################
    # load shower info from inp file
    ########################
    shower = {
            "ID" : showerID,               # shower ID, number of simulation
            "primary" : primarytype,        # primary (electron, pion)
            "energy" : energy,               # EeV
            "zenith" : zen,               # deg (GRAND frame)
            "azimuth" : azim,                # deg (GRAND frame)
            "injection_height" : injh,    # m (injection height in the local coordinate system)
            "task" : task,    # Identification
            "core" : core,    # m, numpy array, core position
            "simulation" : simus # coreas or zhaires
            }
    ####################################
    print("shower", shower)
    
    
    
    
    import astropy
    from astropy.table import Table
    from astropy.table import hstack
    import h5py
    
    i=0
    for ant in glob.glob(path+ending_e):
      while i<1:
        print("\n")
        i+=1
        if simus == 'zhaires':
            ant_number = int(ant.split('/')[-1].split('.trace')[0].split('a')[-1])
            # define path for storage of hdf5 files
            name = path+'/table'+str(ant_number)+'.h5'
            print("Table saved as: ", name)

        if  simus == 'coreas':
            base=os.path.basename(ant)
            # coreas
            ID=(os.path.splitext(base)[0]).split("_")[1] # remove raw 
            ant_number = ID_ant.index(ID)
            # define path for storage of hdf5 files
            name = path+'/../table_'+str(ID)+'.hdf5'
            print("Table saved as: ", name)
            

        ##### read-in output of simulations
        # read in trace from file and store as astropy table, saved as hdf5 file (optional)
        a= load_trace_to_table(path=ant, pos=positions[ant_number], info=shower, content="e", simus=simus, save=name) 
        
        ####### From here on only additional feature and nice-to-know
        ## play around
        #print(a.info)
        #print(a['Ex'])
        #print(a.meta)
        #print(type(a))
        #print(a)

        
        # write astropy table to hdf5 file
        ##a.write(name, path=name, overwrite=True, compression=True, serialize_meta=True) #append=True, 
        #a.write(name, path='efield', format="hdf5",  serialize_meta=True)
        

        
        
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
            
        # example how to add voltage from text file as a column    
        voltage=False
        if voltage:
            ##### Hopefully not needed any more if voltage traces are not stored as txt files in future
            print("Adding voltages")
            
            ### ATTENTION currently electric field added - adjust <ant=path+ending_e>
            print("WARNING: adopt path to voltage trace")
            # read in trace from file and store as astropy table - can be substituted by computevoltage operation
            b= load_trace_to_table(path=ant, pos=positions[ant_number], info=shower, content="v") # read from text files
    
            ## stack electric field and voltage traces
            #from astropy.table import hstack
            #c = hstack([a, b])
            
            # Write to tables in hdf5 file of the efield
            b.write(name, path='voltages', append=True, serialize_meta=True) #append=True -- NOTE: Do I need that
        
        # example VOLATGE COMPUTATION and add to same hdf5 file
        voltage_compute=False
        if voltage_compute:
            from computevoltage import get_voltage, compute_antennaresponse
            from full_chain import run
            efield1=np.array([a['Time'], a['Ex'], a['Ey'], a['Ez']]).T
            # apply only antenna response
            #voltage = compute_antennaresponse(efield1, shower['zenith'], shower['azimuth'], alpha=.0, beta=0. )

            ## apply full chain
            voltage = run(efield1, shower['zenith'], shower['azimuth'], 0, 0, False)
            # load voltage array to table and store in same hdf5 file
            volt_table = _table_voltage(voltage, shower['position'])
            volt_table.write(name, path='voltages', format="hdf5", append=True, serialize_meta=True)

        ######## just testing part and examples how to use astropy tables
        EXAMPLE=False
        if EXAMPLE:
            # read in hdf5 file 
            f=Table.read(name, path="efield")
            g=Table.read(name, path="voltages")
            print(f)
            print(g)
            #print(f)
            #b=f['Ex']*f['Ex']
            #print(b)
            #print(f.meta, f.info)
            #print(f.meta['zenith'], f['Time'].unit)
            #print(g)
            
            ### Just examples how output could handled 
            #summe=f['Ex']+f['Ey']
            #print(summe[-2])
        
      else:
        break
