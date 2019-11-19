################################
#### by A. Zilles
################################

#!/usr/bin/env python

import os
from os.path import split, join, realpath
import numpy as np
import sys
import glob

import logging
logger = logging.getLogger("In_Out")

from astropy import units as u

#===========================================================================================================
#===========================================================================================================
def inputfromtxt(input_file_path):
    # I will move this function as soon as I know to where :)
#===========================================================================================================
    '''
    Get shower parameter from inp file for zhaires simulations
    
    Parameters:
        input_file_path: str
            path of inp file
        
    Returns:
        zen: float
            zenith, deg in GRAND conv.
        azim: float
            azimuth, deg in GRAND conv.
        energy: float
            primary energy in eV
        injh: float
            in meters, 100km for CRs
        primarytype: str
            'proton' or 'iron'
        core: numpy array
            core position in meters
        task: str
            ID of shower
    '''
    
    particule = ['eta','pi+','pi-','pi0','Proton','p','proton','gamma','Gamma','electron','Electron','e-','K+','K-','K0L','K0S','K*+'
    ,'muon+','muon-','Muon+','Muon-','mu+','mu-','tau+','tau-','nu(t)','Positron','positron','e+']


    if os.path.isfile(input_file_path) ==  False:  # File does not exist 
        print('--- ATTENTION: inp-file does not exist')
        #exit()
        
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
                energy = energy
            if unit == "GeV\n":
                energy = energy *1e9
            if unit == "EeV\n":
                energy = energy *1e18
        if 'PrimaryParticle' in line:
            primarytype = str(line.split(' ',-1)[1])
            if primarytype[-1]=='\n':
                primarytype=primarytype[0:-1]
            if primarytype[-1]=='\r':
                primarytype=primarytype[0:-1]
        if 'AddSpecialParticle      RASPASSMulti' in line:
            RASPASSMulti_line = line
        if 'TaskName' in line:
            task = str(line.split(' ',-1)[1])
            if task[-1]=='\n':
                task=task[0:-1]
            if task[-1]=='\r':
                task=task[0:-1]
        if '#Core Position:' in line:
            offset = line.split(' ',-1)
            core = np.array([float(offset[2]), float(offset[3]), float(offset[4])])

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
    try:
        task
    except NameError:
        task = None
    try:
        core
    except NameError:
        core = np.array([0.,0.,0.])


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
    
    # TODO: add observer level -- antenna height to be corrected
    

    #if task:
        #if core.all:
            #return zen,azim,energy,injh,primarytype,core,task
    #if task:
        #return zen,azim,energy,injh,primarytype,task
    #else:
        #return zen,azim,energy,injh,primarytype
        
    return zen,azim,energy,injh,primarytype,core,task

#===========================================================================================================
#===========================================================================================================

#===========================================================================================================
def inputfromtxt_coreas(input_file_path): # still ongoing work
#===========================================================================================================
    '''
    Get shower parameter from inp and reas file for coreas simulations
    
    Parameters:
        input_file_path: str
            path of inp file
        
    Returns:
        zen: float
            zenith, deg in GRAND conv.
        azim: float
            azimuth, deg in GRAND conv.
        energy: float
            primary energy in eV
        injh: not working
            Coreas can not yet accept an injection height
        primarytype: str
            'proton' or 'iron'
        core: numpy array
            core position in meters
        task: str
            ID of shower
    '''
    
    
    if os.path.isfile(input_file_path) ==  False:  # File does not exist 
        print('--- ATTENTION: inp-file does not exist')
        exit()
        
    #datafile = file(input_file_path) # why it is not working...
    datafile = open(input_file_path, 'r') 

    for line in datafile:
        # NOTE: CORSIKA/CoREAS angle = direction of propgation == GRAND conventions
        if 'THETAP' in line:
            zen=float(line.split('    ',-1)[1]) *u.deg# propagation direction
            zen=180*u.deg-zen # to GRAND
        if 'PHIP' in line:
            azim = float(line.split('    ',-1)[1])*u.deg # propagation direction, already in GRAND
        #if 'RASPASSHeight' in line:
            #injh = float(line.split(' ',-1)[2])*u.m
        if 'ERANGE' in line:
            energy = line.split('    ',-1)[1] # GeV by default
            if energy[-1]=='\n':
                energy=energy[0:-1]
            energy = float(energy) *1e9 *u.eV# GeV to eV
        if 'PRMPAR' in line:
            primarytype = str(line.split('    ',-1)[1])
            if primarytype[-1]=='\n':
                primarytype=primarytype[0:-1]
            if primarytype == str(14):
                primarytype ='proton'
            if primarytype == str(5626):
                primarytype ='iron'

    try:
        energy
    except NameError:
        print('No primary energy found in the ZHAireS input text file.')
        exit()
    try:
        primarytype
    except NameError:
        print('ATTENTION: No particle type found')
        primarytype = None
    try:
        injh
    except NameError:
        injh = 100000.e2*u.m #Case of a cosmic for which no injection height is defined in the input file and is then set to 100 km in cm
            
    # Get reas file
    path, reas = os.path.split(input_file_path)
    base = os.path.basename(reas)
    base1 = os.path.splitext(base)
    file_path= path[0:-4]+'/'+base1[0]+".info"
    #file_path= path[0:-4]+'/'+base1[0]+".reas"

    datafile = open(file_path, 'r') 
    for line in datafile:
        if 'TASK' in line:
            task = str(line.split('  ',-1)[1])
            if task[-1]=='\n':
                task=task[0:-1]
            if task[-1]=='\r':
                task=task[0:-1]
            
        if 'CORE' in line:
            #print(line)
            offset = line.split('  ',-1)
            offset[-1]=offset[-1].rstrip()
            core = list([float(offset[1]), float(offset[2]), float(offset[3])])*u.m # in cm to m 
    try:
        task
    except NameError:
        task = None
    try:
        core
    except NameError:
        core = None       
        
        
    if task:
        if core is not None:
            return zen,azim,energy,injh,primarytype,core,task
        else:
            return zen,azim,energy,injh,primarytype,task
    else:
        return zen,azim,energy,injh,primarytype
#===========================================================================================================

def _get_positions_coreas(path):
    '''
    read in antenna positions from Coreas simulations, wrt to sealevel
    
    Parameters:
    datafile: str
        path to folder of run
    
    Returns:
    positions: numpy array
        x,y,z component of antenna positions in meters
    ID_ant: list
        corresponding antenna ID for identification !- [0,1,2,....]
        
        
    NOTE: units assign to positions and slope assuming meters and degree 
    '''
    datafile = open(path, 'r') 
    x_pos1=[]
    y_pos1=[]
    z_pos1=[]
    ID_ant=[]
    #positions=[]
    
    alpha=[]
    beta=[]
    for line in datafile:
    # Coreas
        if 'AntennaPosition =' in line: #list file
            #print(line,line.split('    ',-1) )
            x_pos1.append(float(line.split('  ',-1)[1])/100.) #*u.cm) # cm to m
            y_pos1.append(float(line.split('  ',-1)[2])/100.) #*u.cm) # cm to m
            z_pos1.append(float(line.split('  ',-1)[3])/100.) #*u.cm) # cm to m
            ID_ant.append(str(line.split('  ',-1)[4]))
            alpha.append(0)
            beta.append(0)
        if 'ANTENNA' in line: #info file
            x_pos1.append(float(line.split('  ',-1)[2])) #*u.m) 
            y_pos1.append(float(line.split('  ',-1)[3])) #*u.m) 
            z_pos1.append(float(line.split('  ',-1)[4])) #*u.m) 
            ID_ant.append(str(line.split('  ',-1)[1]))
            alpha.append(float(line.split('  ',-1)[5]))
            beta.append(float(line.split('  ',-1)[6]))
            
            
    x_pos1=np.asarray(x_pos1)
    y_pos1=np.asarray(y_pos1)
    z_pos1=np.asarray(z_pos1)
    positions=np.stack((x_pos1,y_pos1, z_pos1), axis=-1 )*u.m
    slopes=np.stack((alpha, beta), axis=-1 )*u.deg
    #print(ID_ant)    
    
    return positions, ID_ant, slopes

#===========================================================================================================

#def _get_Xmax_coreas(path):


    #for line in glob.glob(path + "*.long"):
    ## Coreas
        #if 'PARAMETERS         =' in line: #list file 
            ##PARAMETERS         =   2.5632E+08 -3.1624E+02  7.2411E+02  3.6003E+01 -9.4986E-03  1.7744E-05
            #Xmax = line[3]
    #return 0
    
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
        
    TODO: read-in hdf5 files and return numpy array
    """

    path = "{:}/a{:}{:}".format(directory, index, suffix)
    with open(path, "r") as f:
        return np.array([list(map(float, line.split())) for line in f])
    
#===========================================================================================================
    
def _table_efield(efield, pos=None, slopes=None, info={}, save=None, ant="/"):
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
    from astropy.table import Table, Column
    
    info.update({'position': pos, 'slopes': slopes})
    
    a = Column(data=efield.T[0],unit=u.ns,name='Time')
    b = Column(data=efield.T[1],unit=u.u*u.V/u.meter,name='Ex')
    c = Column(data=efield.T[2],unit=u.u*u.V/u.meter,name='Ey')
    d = Column(data=efield.T[3],unit=u.u*u.V/u.meter,name='Ez')
    efield_ant = Table(data=(a,b,c,d,), meta=info)
    
    if save:
        efield_ant.write(save, path=ant+'efield', format="hdf5", append=True,  compression=True,serialize_meta=True) #
    #if save is None:
    return efield_ant
    
#===========================================================================================================

def _table_voltage(voltage, pos=None, slopes=None, info={}, save=None, ant="/"):    
    ''' 
    Load voltage trace in table with header info  (numpy array to astropy table)
    
    Parameters
    ---------
    voltage: numpy array
        voltage trace
    pos: numpy array
        position of antenna 
    info: dict
        contains shower info
    
    Returns
    ---------   
    voltage_ant: astropy table
    
    '''
    from astropy.table import Table, Column
    
    info.update({'position': pos, 'slopes': slopes})
    
    a = Column(data=voltage.T[0],unit=u.ns,name='Time')
    b = Column(data=voltage.T[1],unit=u.u*u.V,name='Vx')
    c = Column(data=voltage.T[2],unit=u.u*u.V,name='Vy')
    d = Column(data=voltage.T[3],unit=u.u*u.V,name='Vz')
    voltage_ant = Table(data=(a,b,c,d,), meta=info)
    #print(voltage_ant)
    
    processing_info={'voltage': ('antennaresponse', 'noise', 'filter', 'digitise')}
    if save is not None:
        if 'antennaresponse' in info['voltage']:
            path_tmp=ant+'voltages'
        if 'noise' in info['voltage']:
            path_tmp=ant+'voltages_noise'
        if 'filter' in info['voltage']:
            path_tmp=ant+'filtered'
        if 'digitise' in info['voltage']:
            path_tmp=ant+'voltages_digitise'
        
        voltage_ant.write(save, path=path_tmp, format="hdf5", append=True, compression=True,serialize_meta=True) #
    
    return voltage_ant

#===========================================================================================================

#def load_trace_to_table(path_raw, pos=np.array([0,0,0]), slopes=np.array([0,0]), info=None, content="e", simus="zhaires", save=None, ant="/"):
def load_trace_to_table(path_raw,  pos=None, slopes=None, info={}, content="e", simus="zhaires", save=None, ant="/"):

    """Load data from an electric field trace file to astropy table

   Parameters
   ---------
        path_raw: str 
            path to file -- electric field (.trace) or voltage trace (.dat)
        pos: numpy array, floats
            optional, position of antenna
        info: str
            optional. shower infos
        content: str
            e/efield or v/voltages
        sim: str
            coreas/zhaires, pick the simulation
        save: str 
            optional,path to save a hdf5 file
        

   Returns
   ---------
        astropy table
    """
    
    if content=="efield" or content=="e":
        efield = np.loadtxt(path_raw)
        #zhaires: time in ns and efield in muV/m
        if simus=="coreas": 
            efield.T[0]*=1e9 # s to ns
            ## coreas cgs to SI, V/m to muV/m
            efield.T[1]*=2.99792458e4* 1.e6 
            efield.T[2]*=2.99792458e4* 1.e6 
            efield.T[3]*=2.99792458e4* 1.e6
            
        efield_ant = _table_efield(efield, pos=pos, slopes=slopes, info=info, save=save, ant=ant)
    if content=="voltages" or content=="v":
        voltage = np.loadtxt(path_raw)
        efield_ant = _table_voltage(voltage, pos, slopes=slopes, info=info, save=save, ant=ant)
    

        
    return efield_ant
        
#===========================================================================================================

def _load_eventinfo_fromhdf(path_hdf5):
    """Load data from hdf5 file to numpy array and restore shower info

   Parameters
   ---------
        path_hdf5: str 
            path to hdf5 file     
            
   Returns
   ---------
        shower: list
            shower parameters etc
        position: numpy array in m (cross-check)
            positions of all antennas in array (x,y,z)
        slopes: numpy array in deg
            slopes of all antennas in array (alpha, beta)

    """
    
    from astropy.table import Table

    g=Table.read(path_hdf5, path="/event")
           
    ### Get shower infomation    
    try:        
        ID=g.meta['ID'],               # shower ID, number of simulation
    except:
        ID=None
    try:
        primary=g.meta['primary'],        # primary (electron, pion)
    except:
        ID=None
    try:
        energy=g.meta['energy'],               # EeV
    except:
        energy=None
    try:
        zenith=g.meta['zenith'],               # deg (GRAND frame)
    except:
        zenith=None
    try:
        azimuth=g.meta['azimuth'],                # deg (GRAND frame)
    except:
        azimuth=None
    try:
        injection_height=g.meta['injection_height'],    # m (injection height in the local coordinate system)
    except:
        injection_height=None
    try:
        task=g.meta['task'],    # Identification
    except:
        task=None
    try:
        core=g.meta['core'],    # m, numpy array, core position
    except:
        core=None
    try:
        simulation=g.meta['simulation'] # coreas or zhaires
    except:
        simulation=None

    shower = {
        "ID" : ID[0],               # shower ID, number of simulation
        "primary" : primary[0],        # primary (electron, pion)
        "energy" : energy[0],               # EeV
        "zenith" : zenith[0],               # deg (GRAND frame)
        "azimuth" : azimuth[0],                # deg (GRAND frame)
        "injection_height" : injection_height[0],    # m (injection height in the local coordinate system)
        "task" : task[0],    # Identification
        "core" : core[0],    # m, numpy array, core position
        "simulation" : simulation # coreas or zhaires
        }
    
    try:
        ant_ID=g["ant_ID"]
    except:
        ant_ID=None
    
    try:
        positions=np.array([g["pos_x"], g["pos_y"],g["pos_z"]]).T*g["pos_x"].unit
    except:
        positions=None
    
    try:
        slopes=np.array([g["alpha"], g["beta"]]).T*g["alpha"].unit
    except:
        slopes=None
    
    return shower, ant_ID, positions, slopes

#===========================================================================================================

def _load_path(path_hdf5, path="/analysis"):
    """Load any data from hdf5 file and restores infomation on analysis

   Parameters
   ---------
        path_hdf5: str 
            path to hdf5 file     
        path: str (default set)
            can be any keyword
            
            
   Returns
   ---------
        astropy table 
            returns content of selected path
        info
            returns table infomation
        meta: dict
            
    """
            
    from astropy.table import Table
    
    f=None
    info=None
    meta=None
    try:
        f=Table.read(path_hdf5, path=path)
        try:
            info=f.info
        except:
            print("Info on", path," not availble")
        try:
            meta=f.meta
        except:
            print("Meta on ", path," not availble")
    except:
        logger.warning(path_hdf5, " could not be loaded in _load_path")
    
    return f, info, meta
    
    
#===========================================================================================================

def _load_to_array(path_hdf5, content="efield", ant="/"):
    """Load data from hdf5 file to numpy array and restore traces

   Parameters
   ---------
        path_hdf5: str 
            path to hdf5 file 
        content: str 
            grep efield or voltages traces

   Returns
   ---------
        efield1 or voltage1: numpy array
            containing the electric field or voltage trace trace: time, Ex, Ey, Ez or time, Vx, Vy, Vz
        efield['Time'] or voltage['Time'].unit: str
            unit of time column
        efield['Ex'] or voltage['Vx'].unit: str
            unit of efield or voltage field
        position: numpy array in m (cross-check)
            position of the antenna
        slopes: numpy array in deg
            slopes of the antenna
    """   
    
    from astropy.table import Table
    
    if content=="efield" or content=="e":
        efield=Table.read(path_hdf5, path=ant+"/efield")
        efield1=np.array([efield['Time'], efield['Ex'], efield['Ey'], efield['Ez']]).T
        #print(efield.T[0], efield.T[1])
            
        try:
            position=efield.meta['position']
        except IOError:
            position=None
        try:
            slopes=efield.meta['slopes']
        except IOError:
            slopes=None

        # TODO do we one to return atsropy units...
        return efield1, efield['Time'].unit, efield['Ex'].unit, position, slopes
    
    if content=="voltages" or content=="v":
        try:
            voltage=Table.read(path_hdf5, path=ant+"/voltages")
            voltage1=np.array([voltage['Time'], voltage['Vx'], voltage['Vy'], voltage['Vz']]).T
            #print(voltage.T[0], voltage.T[1])
        except:
            print("Voltages not found")
            voltage=None
            voltage1=None
        
        try:
            position=voltage.meta['position']
        except IOError:
            position=None
        try:
            slopes=voltage.meta['slopes']
        except IOError:
            slopes=None
        
        # TODO do we one to return atsropy units...
        return voltage1, voltage['Time'].unit, voltage['Vx'].unit, position, slopes



#def _load_to_array(path_hdf5, content="efield", ant="/"):
    #"""Load data from hdf5 file to numpy array and restore shower info
       #works only with single antenna files -- TODO modification needed

   #Parameters
   #---------
        #path_hdf5: str 
            #path to hdf5 file 
        #content: str 
            #grep efield or voltages traces

   #Returns
   #---------
        #efield1 or voltage1: numpy array
            #containing the electric field or voltage trace trace: time, Ex, Ey, Ez or time, Vx, Vy, Vz
        #efield['Time'] or voltage['Time'].unit: str
            #unit of time column
        #efield['Ex'] or voltage['Vx'].unit: str
            #unit of efield or voltage field
        #shower: list
            #shower parameters etc
        #position: numpy array in m (cross-check)
        #slopes: numpy array in deg
    #"""   
    #from astropy.table import Table
    
    #if content=="efield" or content=="e":
        #efield=Table.read(path_hdf5, path=ant+"efield")
        #efield1=np.array([efield['Time'], efield['Ex'], efield['Ey'], efield['Ez']])
        ##print(efield.T[0], efield.T[1])
    
        #try:
            #shower = {
                #"ID" : efield.meta['ID'],               # shower ID, number of simulation
                #"primary" : efield.meta['primary'],        # primary (electron, pion)
                #"energy" : efield.meta['energy'],               # EeV
                #"zenith" : efield.meta['zenith'],               # deg (GRAND frame)
                #"azimuth" : efield.meta['azimuth'],                # deg (GRAND frame)
                #"injection_height" : efield.meta['injection_height'],    # m (injection height in the local coordinate system)
                #"task" : efield.meta['task'],    # Identification
                #"core" : efield.meta['core'],    # m, numpy array, core position
                #"simulation" : efield.meta['simulation'] # coreas or zhaires
                #}
        #except:
            #print("NOTE: shower info not contained")
            #shower=None
            
        #try:
            #position=efield.meta['position']
        #except IOError:
            #position=None
        #try:
            #slopes=efield.meta['slopes']
        #except IOError:
            #slopes=None

        ## TODO do we one to return atsropy units...
        #return efield1, efield['Time'].unit, efield['Ex'].unit, shower, position, slopes
    
    #if content=="voltages" or content=="v":
        #try:
            #voltage=Table.read(path_hdf5, path=ant+"voltages")
            #voltage1=np.array([voltage['Time'], voltage['Vx'], voltage['Vy'], voltage['Vz']])
            ##print(voltage.T[0], voltage.T[1])
        #except:
            #print("Voltages not found")
            #voltage=None
            #voltage1=None
            
        #try:
            #shower = {
                #"ID" : voltage.meta['ID'],               # shower ID, number of simulation
                #"primary" : voltage.meta['primary'],        # primary (electron, pion)
                #"energy" : voltage.meta['energy'],               # EeV
                #"zenith" : voltage.meta['zenith'],               # deg (GRAND frame)
                #"azimuth" : voltage.meta['azimuth'],                # deg (GRAND frame)
                #"injection_height" : voltage.meta['injection_height'],    # m (injection height in the local coordinate system)
                #"task" : voltage.meta['task'],    # Identification
                #"core" : voltage.meta['core'],    # m, numpy array, core position
                #"simulation" : voltage.meta['simulation'] # coreas or zhaires
                #}
        #except:
            #print("NOTE: shower info not contained")
            #shower=None
        
        #try:
            #position=voltage.meta['position']
        #except IOError:
            #position=None
        #try:
            #slopes=voltage.meta['slopes']
        #except IOError:
            #slopes=None
        
        
        ## TODO do we one to return atsropy units...
        #return voltage1, voltage['Time'].unit, voltage['Vx'].unit, shower, position, slopes

    

#===========================================================================================================


def load_eventinfo_tohdf(path, showerID, simus, name_all=None):
    """Load data from simulation to hdf5 file

   Parameters
   ---------
        path: str 
            path simulated event set
        showerID: str 
            identifaction string for single event/shower
        simus: str
            coreas/zhaires
        name_all: str (optional)
            path to store the hdf5 file

   Returns
   ---------
        shower: dict
            contains shower parameters and other infos
        ID_ant: numpy array
            antenna ID of whole array
        ID_ant: numpy array
            antenna positions of whole array [x,y,z]
        ID_ant: numpy array
            slopes of whole array [alphha,beta]

        saves hdf5 file with event table and meta data if name_all !=None

   """     
    
    if simus == 'zhaires':
        ####################################### NOTE zhaires --- THOSE HAS TO BE UPDATED
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
            task=None
            core=None 
            
        # correction of shower core
        try:
            positions = positions + np.array([core[0], core[1], 0.])
        except:
            print("positions not corrected for core")
        
        ending_e = "a*.trace"
                

    if  simus == 'coreas':
        #posfile = path +'SIM'+str(showerID)+'.list' # contains not alpha and beta
        posfile = path +'SIM'+str(showerID)+'.info' # contains original ant ID , positions , alpha and beta
        positions, ID_ant, slopes = _get_positions_coreas(posfile)
        #print(positions, ID_ant, slopes)
            
        inputfile = path+'/inp/SIM'+showerID+'.inp'
        zen,azim,energy,injh,primarytype,core,task = inputfromtxt_coreas(inputfile)
           
        # correction of shower core
        try:
            positions = positions + np.array([core[0], core[1], 0.*u.m])
        except:
            logger.debug("No core position information availble")
    
        #----------------------------------------------------------------------   

        
    # load shower info from inp file via dictionary
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
    #print("shower", shower)
    logger.info("Shower summary: " + str(shower))
        

    
    if name_all is not None:
        from astropy.table import Table, Column
        a1 = Column(data=np.array(ID_ant), name='ant_ID')
        b1 = Column(data=positions.T[0], unit=u.m, name='pos_x')
        c1 = Column(data=positions.T[1], unit=u.m, name='pos_y')
        d1 = Column(data=positions.T[2], unit=u.m, name='pos_z')  #u.eV, u.deg
        e1 = Column(data=slopes.T[0], unit=u.deg, name='alpha')
        f1 = Column(data=slopes.T[1], unit=u.deg, name='beta') 
        event_info = Table(data=(a1,b1,c1,d1,e1,f1,), meta=shower) 
        event_info.write(name_all, path='event', format="hdf5", append=True,  compression=True, serialize_meta=True)
        print("Event info saved in: ", name_all)

    return shower, ID_ant, positions, slopes


