"""
My first try to write a python class

Shower class: shall contain all infos on the shower
"""

import numpy as np

class shower:
    ''' info on shower parameter 
        
        showerID: str
        primary: str
        energy: float in eV
        zenith: float in deg, GRAND
        azmiuth: float in deg, GRAND
        
        simulations: str
        
        recoXmax: float in g/cm2
        recoenergy: float in eV
        recozenith: float in deg, GRAND
        recoazimuth: float in deg, GRAND
        
    '''
    def __init__(self):
        self.max_length = 1
        
        self.showerID = [] # attribute references
        self.primary = []
        self.energy = []
        self.zenith = []    # How to limit to one number per list .....
        self.azimuth = [] #instantiation
        self.injectionheight = []
        
        # simulated shower
        self.simulation = []
        
        # reconstructed shower
        self.recoXmax = []
        self.recoenergy = []
        self.recozenith = []
        self.recoazimuth = []
    
    def add_showerID(self, showerID):
        if len(self.showerID)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add showerID=", str(showerID)," -- already set to: ", self.showerID[0])
        elif showerID is None:
            print("Warning: Primary not given")
        else:
            self.showerID.append(showerID)
        
    def get_showerID(self):
        if len(self.showerID)==0:
            print("Warning: no showerID set")
        else:
            return self.showerID[0] # str
        
    
    def add_primary(self, primary):
        if len(self.primary)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add primary=", str(primary)," -- already set to: ", self.primary[0])
        elif primary is None:
            print("Warning: Primary not given")
        else:
            self.primary.append(primary)
        
    def get_primary(self):
        if len(self.primary)==0:
            print("Warning for shower ", str(self.showerID),": No primary set")
        else:
            return self.primary[0] # str
        
    def add_energy(self, energy):
        if len(self.energy)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add energy=", str(energy)," -- already set to: ", self.energy[0])
        elif energy is None:
            print("Warning for shower ", str(self.showerID),": Energy not given")
        else:
            self.energy.append(energy)
        
    def get_energy(self):
        if len(self.energy)==0:
            print("Warning for shower ", str(self.showerID),": No energy set")
        else:
            return self.energy[0]  
        
    def add_zenith(self, zenith):
        if len(self.zenith)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add zenith=", str(zenith)," -- already set to: ", self.zenith[0])
        elif zenith is None:
            print("Warning for shower ", str(self.showerID),": Zenith not given")
        else:
            self.zenith.append(zenith)

    def get_zenith(self):
        if len(self.zenith)==0:
            print("Warning for shower ", str(self.showerID),": No zenith set")
        else:
            return self.zenith[0]
        
    def add_azimuth(self, azimuth):
        if len(self.azimuth)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add azimuth=", str(azimuth)," -- already set to: ", self.azimuth[0])
        elif azimuth is None:
            print("Warning for shower ", str(self.showerID),": Azimuth not given")
        else:
            self.azimuth.append(azimuth)
        
    def get_azimuth(self):
        if len(self.azimuth)==0:
            print("Warning for shower ", str(self.showerID),": No azimuth set")
        else:
            return self.azimuth[0]
        
    def add_injectionheight(self, injectionheight):
        if len(self.injectionheight)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add injectionheight=", str(injectionheight)," -- already set to: ", self.injectionheight[0])
        elif injectionheight is None:
            print("Warning for shower ", str(self.showerID),": Azimuth not given")
        else:
            self.injectionheight.append(injectionheight)
        
    def get_injectionheight(self):
        if len(self.injectionheight)==0:
            print("Warning for shower ", str(self.showerID),": No injectionheight set")
        else:
            return self.injectionheight[0]
        
        
    def add_all(self,showerID, primary, energy, zenith, azimuth, injectionheight):
        self.add_showerID(showerID)
        self.add_primary(primary)
        self.add_energy(energy)
        self.add_zenith(zenith)
        self.add_azimuth(azimuth)
        self.add_injectionheight(injectionheight)
        
    def get_all(self):
        print(self.get_showerID(), self.get_primary(), self.get_energy(), self.get_zenith(), self.get_azimuth(), self.get_injectionheight())
        
        
    #def get_Xmax(self):
        #computes Xmax for shower
        
    
    #def get_direction(self):
    #np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])
 
#=================================================        
 
### simulated shower 
class sim_shower(shower):
    ''' info on simulations '''
    def add_simulation(self, simulation):
        self.simulation.append(simulation)
        
    def get_simulation(self):
        print(self.simulation)
        
#=================================================        

### reconstructed event
class reco_shower(shower):
    ''' info on reconstructed values '''
    # missing reco injectionheight
   
    def add_recoenergy(self, recoenergy):
        if len(self.recoenergy)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add recoenergy=", str(recoenergy)," -- already set to: ", self.recoenergy[0])
        elif recoenergy is None:
            print("Warning for shower ", str(self.showerID),": recoenergy not given")
        else:
            self.recoenergy.append(recoenergy)
        
    def get_recoenergy(self):
        if len(self.recoenergy)==0:
            print("Warning for shower ", str(self.showerID),": No recoenergy set")
        else:
            return self.recoenergy[0]  
        
    def add_recoXmax(self, recoXmax):
        if len(self.recoXmax)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add recoXmax=", str(recoXmax)," -- already set to: ", self.recoXmax[0])
        elif recoXmax is None:
            print("Warning for shower ", str(self.showerID),": recoXmax not given")
        else:
            self.recoXmax.append(recoXmax)
        
    def get_recoXmax(self):
        if len(self.recoXmax)==0:
            print("Warning for shower ", str(self.showerID),": No recoXmax set")
        else:
            return self.recoXmax[0]  
        
    def add_recozenith(self, recozenith):
        if len(self.recozenith)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add recozenith=", str(recozenith)," -- already set to: ", self.recozenith[0])
        elif recozenith is None:
            print("Warning for shower ", str(self.showerID),": recozenith not given")
        else:
            self.recozenith.append(recozenith)
        
    def get_recozenith(self):
        if len(self.recozenith)==0:
            print("Warning for shower ", str(self.showerID),": No recozenith set")
        else:
            return self.recozenith[0]  
        
    def add_recoazimuth(self, recoazimuth):
        if len(self.recoazimuth)==self.max_length:
            print("Warning for shower ", str(self.showerID),": Can't add recoazimuth=", str(recoazimuth)," -- already set to: ", self.recoazimuth[0])
        elif recoazimuth is None:
            print("Warning for shower ", str(self.showerID),": recoazimuth not given")
        else:
            self.recoazimuth.append(recoazimuth)
        
    def get_recoazimuth(self):
        if len(self.recoazimuth)==0:
            print("Warning for shower ", str(self.showerID),": No recoazimuth set")
        else:
            return self.recoazimuth[0]  
       
       
    def add_recoall(self, recoXmax, energy, zenith, azimuth):
        self.add_recoXmax(recoXmax)
        self.add_recoenergy(energy)
        self.add_recozenith(zenith)
        self.add_recoazimuth(azimuth)
        
    def get_recoall(self):
        print(self.get_recoXmax(), self.get_recoenergy(), self.get_recozenith(), self.get_recoazimuth())
        
        
#=================================================        


def loadInfo_toShower(name, info=None):
    #load meta info from hdf5 file to class
    try:
        showerID = info["ID"]
    except IOError:
        showerID = None
    
    try:
        primary = info["primary"]
    except IOError:
        primary = None
        
    try:
        energy = info["energy"]
    except IOError:
        energy = None
    
    try:
        zenith = info["zenith"]
    except IOError:
        zenith = None
    
    try:
        azimuth = info["azimuth"]
    except IOError:
        azimuth = None
    
    try:
        injectionheight = info["injection_height"]
    except IOError:
        injectionheight = None

    try:
        simulation = info["simulation"]
    except IOError:
        simulation = None
    
    name.add_all(showerID, primary, energy, zenith, azimuth, injectionheight)
    name.add_simulation(simulation)




#############
## TESTING ##
#############
    
#print("\n shower 1")
#sh1 = shower()
#sh1.add_showerID("1")
#sh1.get_zenith()
#print(type(sh1.get_showerID()), sh1.get_showerID())

#sh1.add_all("1", 'electron', 1e17, 93., 187, 2800)
#print(sh1.get_zenith())
#sh1.add_zenith(65)

#sh1.add_azimuth(180)
#print(sh1.get_azimuth())

#print("\n shower 2")
#sh2 = sim_shower()
#sh2.add_all("2", "pion", 3e17, 85., 10, None)
#sh2.add_simulation("coreas")
#sh2.get_simulation()
#sh2.get_all()

#print("\n shower 3")
#sh3 = reco_shower()
#sh3.add_recoall(865, 5e19, 85., 100, )
#sh3.get_recoall()

#import os
#from os.path import split, join, realpath
#import sys
#from in_out import inputfromtxt, _get_positions_coreas, inputfromtxt_coreas, load_trace_to_table

#from astropy.table import Table

#f = Table.read('/home/laval1NS/zilles/CoREAS/GP300-v3/000002/table_a92.hdf5', path="efield") #txt.T[0]:time in ns, txt.T[1]: North-South, txt.T[2]: East-West, txt.T[3]: Up , all electric field in muV/m
#trace = np.array([f['Time'], f['Ex'], f['Ey'], f['Ez']]).T
#unit='muV/m'

#print(f.meta)
#testshower = sim_shower()
#loadInfo_toShower(testshower, f.meta)
#testshower.get_all()

#g = Table.read('/home/laval1NS/zilles/CoREAS/GP300-v3/000002/table_a236.hdf5', path="efield")

#shbla = sim_shower()
#loadInfo_toShower(shbla, g.meta)

