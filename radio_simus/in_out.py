################################
#### by A. Zilles
################################

#!/usr/bin/env python

import numpy as np
import sys
import glob
from astropy.table import Table
from astropy.table import hstack

#===========================================================================================================
#===========================================================================================================
def inputfromtxt(input_file_path):
    # I will move this function as soon as I know to where :)
#===========================================================================================================
    particule = ['eta','pi+','pi-','pi0','Proton','p','proton','gamma','Gamma','electron','Electron','e-','K+','K-','K0L','K0S','K*+'
    ,'muon+','muon-','Muon+','Muon-','mu+','mu-','tau+','tau-','nu(t)','Positron','positron','e+']

    #datafile = file(input_file_path) # why it is not working...
    datafile = open(input_file_path, 'r') 

    for line in datafile:
        if 'PrimaryZenAngle' in line:
            zen=float(line.split(' ',-1)[1])
            zen = 180-zen  #conversion to GRAND convention i.e. pointing towards antenna/propagtion direction
        if 'PrimaryAzimAngle' in line:
            azim = float(line.split(' ',-1)[1])+180 #conversion to GRAND convention i.e. pointing towards antenna/propagtion direction
            if azim>=360:
                azim= azim-360
        if 'RASPASSHeight' in line:
            injh = float(line.split(' ',-1)[2])
        if 'PrimaryEnergy' in line:
            energy = float(line.split(' ',-1)[1])
            unit= str(line.split(' ',-1)[2])
            if unit == "eV\n":
                energy = energy *1e-18
            if unit == "GeV\n":
                energy = energy *1e-9
        if 'PrimaryParticle' in line:
            primarytype = str(line.split(' ',-1)[1])
            if primarytype[-1]=='\n':
                primarytype=primarytype[0:-1]
        if 'AddSpecialParticle      RASPASSMulti' in line:
            RASPASSMulti_line = line

    try:
        injh
    except NameError:
        injh = 100000. #Case of a cosmic for which no injection height is defined in the input file and is then set to 100 km by ZHAireS
    try:
        energy
    except NameError:
        print('No primary energy found in the ZHAireS input text file.')
        exit()
    try:
        primarytype
    except NameError:
        primarytype = None


    if primarytype=='RASPASSMulti':
        tmp = RASPASSMulti_line.split(' ',-1)
        if tmp[-1][-1]=='\n':
            tmp[-1]=tmp[-1][0:-1]
        prod = [x for x in particule if x in set(tmp)]
        ind_prod = np.array([tmp.index(x) for x in prod],dtype=int)
        Wprod = [float(tmp[ind]) for ind in ind_prod+1]
        primarytype = prod[np.argmax(Wprod)]

    if primarytype=='Proton' or primarytype=='K+' or primarytype=='K-' or primarytype=='K0L' or primarytype=='K0S' or primarytype=='K*+':
        primarytype='proton'
    elif primarytype=='gamma' or primarytype=='Gamma' or primarytype=='Electron':
        primarytype='electron'
    elif primarytype=='pi0' or primarytype=='pi-' or primarytype=='pi+':
        primarytype='pion'

    return zen,azim,energy,injh,primarytype


#===========================================================================================================
#===========================================================================================================

def load_trace(directory, index, suffix=".trace"):
    """Load data from a trace file

   Parameters
   ---------
        directory: str 
            path to file
        index: ind
            index number of antenna
        suffix: str 
            optional, suffix of file

   Returns
   ---------
        numpy array
    """

    path = "{:}/a{:}{:}".format(directory, index, suffix)
    with open(path, "r") as f:
        return numpy.array([list(map(float, line.split())) for line in f])
    
#===========================================================================================================
    
def _table_efield(efield, pos, info={}):
    ''' 
    Load electric field trace in table with header info (numpy array to astropy table)
    
    Parameters
    ---------
    efield: numpy array
        electric field trace
    pos: numpy array
        position of antenna 
    info: dict
        contains shower info
    
    Returns
    ---------   
    efield_ant: astropy table
        
    '''
    
    info.update({'position': pos})
    efield_ant = Table(efield, names=('Time', 'Ex', 'Ey', 'Ez'), meta=info)
    efield_ant['Time'].unit= 'ns'
    efield_ant['Ex'].unit= 'muV/m'
    efield_ant['Ey'].unit= 'muV/m'
    efield_ant['Ez'].unit= 'muV/m'
    return efield_ant
    
#===========================================================================================================

def _table_voltage(voltage, pos, info={}):    
    ''' 
    Load voltage trace in table with header info  (numpy array to astropy table)
    
    Parameters
    ---------
    efield: numpy array
        voltage trace
    pos: numpy array
        position of antenna 
    info: dict
        contains shower info
    
    Returns
    ---------   
    voltage_ant: astropy table
    
    '''
    info.update({'position': pos})
    voltage_ant = Table(voltage, names=('Time', 'Vx', 'Vy', 'Vz'), meta=info)
    voltage_ant['Time'].unit= 's'
    voltage_ant['Vx'].unit= 'muV'
    voltage_ant['Vy'].unit= 'muV'
    voltage_ant['Vz'].unit= 'muV'
    return voltage_ant

#===========================================================================================================

def load_trace_to_table(directory, index, pos=np.array([0,0,0]), info=None, suffix=".trace"):
    """Load data from an electric field trace file to astropy table

   Parameters
   ---------
        directory: str 
            path to file -- electric field (.trace) or voltage trace (.dat)
        index: ind
            index number of antenna
        pos: numpy array, floats
            optional, position of antenna
        suffix: str 
            optional, suffix of file

   Returns
   ---------
        astropy table
    """
    
    if suffix==".trace":
        path = "{:}/a{:}{:}".format(directory, index, suffix)
        efield = np.loadtxt(path)
        efield_ant = _table_efield(efield, pos, info)
    if suffix==".dat":
        path = "{:}/out_{:}{:}".format(directory, index, suffix)
        voltage = np.loadtxt(path)
        efield_ant = _table_voltage(voltage, pos)

        
    return efield_ant
        

    
#===========================================================================================================
#===========================================================================================================

if __name__ == '__main__':

    if ( len(sys.argv)<1 ):
        print("""
        Example how to read in an electric field or voltage trace, load it to a table, write it into a hdf5 file and read the file.
        -- produces hdf5 files
        
        Usage for full antenna array:
            python full_chain.py [path to folder]
        example: python in_out.py ./
            
        """)
        sys.exit(0)


    # path to folder containing traces
    path = sys.argv[1]
    
    # Get the antenna positions from file
    positions = np.loadtxt(path+"antpos.dat")
        
    # Get shower info
    showerID=str(path).split("/")[-2] # should be equivalent to folder name
    inputfile = path+'/inp/'+showerID+'.inp'
    
    zen,azim,energy,injh,primarytype = inputfromtxt(inputfile)
    ########################
    # TODO: load shower info from inp file
    ########################
    shower = {
        "ID" : showerID,               # shower ID, number of simulation
        "primary" : primarytype,        # primary (electron, pion)
        "energy" : energy,               # EeV
        "zenith" : zen,               # deg (GRAND frame)
        "azimuth" : azim,                # deg (GRAND frame)
        "injection_height" : injh    # m (injection height in the local coordinate system) 
        }
    #print(shower)
    
    for ant in glob.glob(path+'*.trace'):

        ant_number = int(ant.split('/')[-1].split('.trace')[0].split('a')[-1])
        
        ##### read-in output of simulations
        print("Getting electric field traces")
        
        # read in trace from file and store as astropy table
        a= load_trace_to_table(path, ant_number, pos=positions[ant_number], info=shower, suffix=".trace")
        #print(a.info)
        #print(a['Ex'])
        #print(a.meta)
        
        # define a path, NOTE: I havent fully understood the path thingy
        name = path+'/table'+str(ant_number)+'.hdf5'
        
        # write astropy table to hdf5 file
        a.write(name, path=name, overwrite=True, serialize_meta=True) #append=True, 
        
        # read in hdf5 file 
        read_a = Table.read(name, path=name)
        #print(read_a.meta['zenith'])
        
        ### Just examples how output could handled 
        #summe=a['Ex']+a['Ey']
        #print(summe[-2])
        
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
            

        ##### Hopefully not needed any more if voltage traces are not stored as txt files in future
        print("Adding voltages")
        
        # read in trace from file and store as astropy table - can be substituted by computevoltage operation
        b= load_trace_to_table(path, ant_number, pos=positions[ant_number], info=shower, suffix=".dat")
        
        # stack electric field and voltage traces
        c = hstack([a, b])
        
        # Write to tables in hdf5 file
        c.write(name, path=name, overwrite=True, serialize_meta=True) #append=True -- NOTE: Do I need that

        # read in hdf5 file 
        read_c = Table.read(name, path=name)
        #print(read_b.meta, read_b.info)
        print(read_c)
