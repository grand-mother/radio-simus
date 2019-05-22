################################
#### by A. Zilles, last update: 22 May 2019
################################


import numpy as np
import sys
import glob
from astropy.table import Table
from astropy.table import hstack

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
    voltage_ant['Time'].unit= 'ns'
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
    
    ########################
    # TODO: load shower info from inp file
    ########################
    shower = {
        "ID" : "000001",               # shower ID, number of simulation
        "primary" : "electron",        # primary (electron, pion)
        "energy" : 0.96,               # EeV
        "zenith" : 89.5,               # deg (GRAND frame)
        "azimuth" : 0.,                # deg (GRAND frame)
        "injection_height" : 2000.,    # m (injection height in the local coordinate system)
        "altitude" : 2000. }  
    
    for ant in glob.glob(path+'*.trace'):

        ant_number = int(ant.split('/')[-1].split('.trace')[0].split('a')[-1])
        
        ##### read-in output of simulations
        print("Getting electric field traces")
        
        # read in trace from file and store as astropy table
        a= load_trace_to_table(path, ant_number, pos=positions[ant_number], info=shower, suffix=".trace")
        print(a.info)
        #print(a['Ex'])
        #print(a.meta)
        # TODO: How do I get the raw/pure data without the header info for a further analysis....
        
        # define a path, NOTE: I havent fully understood the path thingy
        name = path+'/table'+str(ant_number)+'.hdf5'
        
        # write astropy table to hdf5 file
        a.write(name, path=name, overwrite=True, serialize_meta=True) #append=True, 
        
        # read in hdf5 file 
        read_a = Table.read(name, path=name)
        
        
        
        
        ##### Hopefully not needed any more if voltage traces are not stored as txt files in future
        print("Adding voltages")
        
        # read in trace from file and store as astropy table - can be substituted by computevoltage operation
        b= load_trace_to_table(path, ant_number, pos=positions[ant_number], info=shower, suffix=".out")
        
        # stack electric field and voltage traces
        c = hstack([a, b])
        
        # Write to tables in hdf5 file
        b.write(name, path=name, overwrite=True, serialize_meta=True) #append=True -- NOTE: Do I need that

        # read in hdf5 file 
        read_b = Table.read(name, path=name)
        #print(read_b.meta, read_b.info)
        
