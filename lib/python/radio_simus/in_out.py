################################
#### by A. Zilles
################################

#!/usr/bin/env python

import os
from os.path import split, join, realpath
import numpy as np
import sys
import glob

from astropy import units as u
#muV_m = u.u * u.V / u.m

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
            zen=float(line.split('    ',-1)[1]) # propagation direction
            zen=180-zen # to GRAND
        if 'PHIP' in line:
            azim = float(line.split('    ',-1)[1]) # propagation direction, already in GRAND
        #if 'RASPASSHeight' in line:
            #injh = float(line.split(' ',-1)[2])
        if 'ERANGE' in line:
            energy = line.split('    ',-1)[1] # GeV by default
            if energy[-1]=='\n':
                energy=energy[0:-1]
            energy = float(energy) *1e9 # GeV to eV
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
        injh = 100000.e2 #Case of a cosmic for which no injection height is defined in the input file and is then set to 100 km in cm
            
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
            if offset[-1]=='\n':
                offset=offset[0:-1]
            core = np.array([float(offset[1]), float(offset[2]), float(offset[3])]) # in cm to m 
    try:
        task
    except NameError:
        task = None
    try:
        core
    except NameError:
        core = None       
        
        
    if task:
        if core.all()!=None:
            return zen,azim,energy,injh,primarytype,core,task
    if task:
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
    positions=np.stack((x_pos1,y_pos1, z_pos1), axis=-1 )
    slopes=np.stack((alpha, beta), axis=-1 )    
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
    """

    path = "{:}/a{:}{:}".format(directory, index, suffix)
    with open(path, "r") as f:
        return numpy.array([list(map(float, line.split())) for line in f])
    
#===========================================================================================================
    
def _table_efield(efield, pos, slopes=None, info={}):
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
    #efield_ant = Table(efield, names=('Time', 'Ex', 'Ey', 'Ez'), meta=info)
    #efield_ant['Time'].unit= 'ns'
    #efield_ant['Ex'].unit= 'muV/m'
    #efield_ant['Ey'].unit= 'muV/m'
    #efield_ant['Ez'].unit= 'muV/m'
    
    a = Column(data=efield.T[0],unit=u.ns,name='Time')
    b = Column(data=efield.T[1],unit=u.u*u.V/u.meter,name='Ex')
    c = Column(data=efield.T[2],unit=u.u*u.V/u.meter,name='Ey')
    d = Column(data=efield.T[3],unit=u.u*u.V/u.meter,name='Ez')
    efield_ant = Table(data=(a,b,c,d,), meta=info)
    
    return efield_ant
    
#===========================================================================================================

def _table_voltage(voltage, pos=None, slopes=None, info={}):    
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
    
    #print(voltage.T[0])
    
    info.update({'position': pos, 'slopes': slopes})
    #voltage_ant = Table(voltage, names=('Time', 'Vx', 'Vy', 'Vz'), meta=info)
    #voltage_ant['Time'].unit= 's'
    #voltage_ant['Vx'].unit= 'muV'
    #voltage_ant['Vy'].unit= 'muV'
    #voltage_ant['Vz'].unit= 'muV'
    
    a = Column(data=voltage.T[0],unit=u.ns,name='Time')
    b = Column(data=voltage.T[1],unit=u.u*u.V,name='Vx')
    c = Column(data=voltage.T[2],unit=u.u*u.V,name='Vy')
    d = Column(data=voltage.T[3],unit=u.u*u.V,name='Vz')
    voltage_ant = Table(data=(a,b,c,d,), meta=info)
    #print(voltage_ant)
    
    return voltage_ant

#===========================================================================================================

def load_trace_to_table(path, pos=np.array([0,0,0]), slopes=np.array([0,0]), info=None, content="e", simus="zhaires", save=None, ant="/"):

    """Load data from an electric field trace file to astropy table

   Parameters
   ---------
        path: str 
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
    
    #if suffix==".trace":
    if content=="efield" or content=="e":
        #path = "{:}/a{:}{:}".format(directory, index, suffix)
        efield = np.loadtxt(path)
        #zhaires: time in ns and efield in muV/m
        if simus=="coreas": 
            efield.T[0]*=1e9 # s to ns
            ## coreas cgs to SI, V/m to muV/m
            efield.T[1]*=2.99792458e4* 1.e6 
            efield.T[2]*=2.99792458e4* 1.e6 
            efield.T[3]*=2.99792458e4* 1.e6
            
        efield_ant = _table_efield(efield, pos=pos, slopes=slopes, info=info)
    #if suffix==".dat":
    if content=="voltages" or content=="v":
        #path = "{:}/out_{:}{:}".format(directory, index, suffix)
        voltage = np.loadtxt(path)
        efield_ant = _table_voltage(voltage, pos)
    
    if save:
        if content=="efield" or content=="e":
            efield_ant.write(save, path=ant+'efield', format="hdf5", append=True,  compression=True,serialize_meta=True) #
        if content=="voltages" or content=="v":
            efield_ant.write(save, path=ant+'voltages', format="hdf5", append=True, compression=True,serialize_meta=True) #
        
    return efield_ant
        

def _load_to_array(path_hdf5, content="efield"):
    """Load data from hdf5 file to numpy array and restore shower info

   Parameters
   ---------
        path_hdf5: str 
            path to hdf5 file 
        content: str 
            grep efield or voltages traces

   Returns
   ---------
        voltage1: numpy array
            containing the electric field trace: time, Ex, Ey, Ez
        voltage['Time'].unit
            unit of time column
        voltage['Ex'].unit
            unit of electric field
        shower: list
            shower parameters etc
        position: numpy array
        slopes: numpy array
    """   
    from astropy.table import Table
    
    if content=="efield" or content=="e":
        efield=Table.read(path_hdf5, path="efield")
        efield1=np.array([efield['Time'], efield['Ex'], efield['Ey'], efield['Ez']])
        #print(efield.T[0], efield.T[1])
    
        try:
            shower = {
                "ID" : efield.meta['ID'],               # shower ID, number of simulation
                "primary" : efield.meta['primary'],        # primary (electron, pion)
                "energy" : efield.meta['energy'],               # EeV
                "zenith" : efield.meta['zenith'],               # deg (GRAND frame)
                "azimuth" : efield.meta['azimuth'],                # deg (GRAND frame)
                "injection_height" : efield.meta['injection_height'],    # m (injection height in the local coordinate system)
                "task" : efield.meta['task'],    # Identification
                "core" : efield.meta['core'],    # m, numpy array, core position
                "simulation" : efield.meta['simulation'] # coreas or zhaires
                }
        except:
            print("NOTE: shower info not contained")
            shower=None
            
        try:
            position=efield.meta['position']
        except IOError:
            position=None
        try:
            slopes=efield.meta['slopes']
        except IOError:
            slopes=None

        return efield1, efield['Time'].unit, efield['Ex'].unit, shower, position, slopes
    
    if content=="voltages" or content=="v":
        try:
            voltage=Table.read(path_hdf5, path="voltages")
            voltage1=np.array([voltage['Time'], voltage['Vx'], voltage['Vy'], voltage['Vz']])
            #print(voltage.T[0], voltage.T[1])
        except:
            print("Voltages not found")
            voltage=None
            voltage1=None
            
        try:
            shower = {
                "ID" : voltage.meta['ID'],               # shower ID, number of simulation
                "primary" : voltage.meta['primary'],        # primary (electron, pion)
                "energy" : voltage.meta['energy'],               # EeV
                "zenith" : voltage.meta['zenith'],               # deg (GRAND frame)
                "azimuth" : voltage.meta['azimuth'],                # deg (GRAND frame)
                "injection_height" : voltage.meta['injection_height'],    # m (injection height in the local coordinate system)
                "task" : voltage.meta['task'],    # Identification
                "core" : voltage.meta['core'],    # m, numpy array, core position
                "simulation" : voltage.meta['simulation'] # coreas or zhaires
                }
        except:
            print("NOTE: shower info not contained")
            shower=None
        
        try:
            position=voltage.meta['position']
        except IOError:
            position=None
        try:
            slopes=voltage.meta['slopes']
        except IOError:
            slopes=None
        
        return voltage1, voltage['Time'].unit, voltage['Ex'].unit, shower, position, slopes

    


