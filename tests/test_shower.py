# -*- coding: utf-8 -*-
"""
Unit tests for the radio_simus.shower module

Usage: python3.7 container_examples/test_shower.py


To be moved in tests folder
"""

import unittest
import sys
import numpy as np

import astropy.units as u

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
#from radio_simus.shower import * 
from radio_simus.shower import Shower , SimulatedShower, ReconstructedShower
##from framework import git ## TODO: not yet linked


class ShowerTest(unittest.TestCase):
    """Unit tests for the version module"""
    
    #def _init_(self): # Constructor
    shower1 = Shower()


    def test_showerID(self):
        print("Start _showerID test")
        ID = []
        for i in range(1): #(2):
            shower_ID=None
            ID.append("name"+str(i))            
            # add ID to shower ---  can only get one ID per shower
            self.shower1.showerID = str("name"+str(i))    
            # get ID from shower
            shower_ID = self.shower1.showerID
                        
            if i==0:
                self.assertEqual(ID[0], shower_ID) #shall be equal
                self.assertEqual(len(ID[0]), len(shower_ID))
            #if i==1: # test whether variable can be overwritten
                #self.assertNotEqual(ID[1], shower_ID) # shall be not equal
                #self.assertEqual(len(ID[1]), len(shower_ID))
                
        print("Finish _showerID test")
        

    def test_primary(self):
        print("Start _primary test")
        prim = []

        prim.append("proton")            
        # add primary to shower ---  can only get one primary element per shower
        self.shower1.primary = prim[-1]    
        # get primary from shower
        primary = self.shower1.primary
        self.assertEqual(prim[0], primary) #shall be equal
        
        ## test whether variable can be overwritten      
        #prim.append("electron")            
        ## add primary to shower ---  can only get one primary element per shower
        #self.shower1.add_primary(prim[-1])    
        ## get primary from shower
        #primary = self.shower1.get_primary()        
        #self.assertNotEqual(prim[1], primary) # shall be not equal
                
        print("Finish _primary test")
        
        
    def test_energy(self):
        print("Start _energy test")
        ener = []

        ener.append(1e18* u.eV )            
        # add energy to shower ---  can only get one primary element per shower
        self.shower1.energy = ener[-1] 
        # get energy from shower
        energy = self.shower1.energy
        self.assertEqual(ener[0], energy) #shall be equal       
        
        ## test whether variable can be overwritten              
        #ener.append(3.12e17)            
        ## add energy to shower ---  can only get one primary element per shower
        #self.shower1.add_energy(ener[-1])    
        ## get energy from shower
        #energy = self.shower1.get_energy()        
        #self.assertNotEqual(ener[1], energy) # shall be not equal
                
        print("Finish _energy test")
        
        
    def test_zenith(self):
        print("Start _zenith test")
        zen = []

        zen.append(70.* u.deg)            
        # add zenithto shower ---  can only get one primary element per shower
        self.shower1.zenith = zen[-1]
        # get zenith from shower
        zenith = self.shower1.zenith
        self.assertEqual(zen[0], zenith) #shall be equal       

        ## test whether variable can be overwritten         
        #zen.append(83.)            
        ## add energy to shower ---  can only get one primary element per shower
        #self.shower1.add_zenith(zen[-1])    
        ## get energy from shower
        #zenith = self.shower1.get_zenith()        
        #self.assertNotEqual(zen[1], zenith) # shall be not equal
                
        print("Finish _zenith test")
                

    def test_azimuth(self):
        print("Start _azimuth test")
        azim = []

        azim.append(100.* u.deg)            
        # add azimuthto shower ---  can only get one primary element per shower
        self.shower1.azimuth = azim[-1]    
        # get azimuth from shower
        azimuth = self.shower1.azimuth
        self.assertEqual(azim[0], azimuth) #shall be equal       
        
        ## test whether variable can be overwritten         
        #azim.append(185.)            
        ## add energy to shower ---  can only get one primary element per shower
        #self.shower1.add_azimuth(azim[-1])    
        ## get energy from shower
        #azimuth = self.shower1.get_azimuth()        
        #self.assertNotEqual(azim[1], azimuth) # shall be not equal
                
        print("Finish _azimuth test")        
        
        
    def test_injectionheight(self):
        print("Start _injectionheight test")
        injh = []

        injh.append(1500.* u.m)            
        # add injectionheightto shower ---  can only get one primary element per shower
        self.shower1.injectionheight = injh[-1]    
        # get injectionheight from shower
        injectionheight = self.shower1.injectionheight
        self.assertEqual(injh[0], injectionheight) #shall be equal       
        
        #injh.append(2000.)            
        ## add energy to shower ---  can only get one primary element per shower
        #self.shower1.add_injectionheight(injh[-1])    
        ## get energy from shower
        #injectionheight = self.shower1.get_injectionheight()        
        #self.assertNotEqual(injh[1], injectionheight) # shall be not equal
                
        print("Finish _injectionheight test")             
        
        
    def test_trigger(self):
        print("Start _trigger test")
        tr = []

        tr.append((1, "yes", 45* u.mV))            
        # add injectionheightto shower ---  can only get one primary element per shower
        self.shower1.trigger = tr[-1]    
        # get trigger from shower
        trigger = self.shower1.trigger
        self.assertEqual(tr[0], trigger) #shall be equal       
        
        print("Finish _trigger test")    



    shower2 = Shower(azimuth = 60.* u.deg, zenith = 45.* u.deg)

    def test_direction(self):
        print("Start _direction test")
        #shower2 = Shower(azimuth = 60.* u.deg, zenith = 45.* u.deg)
        
        azi = self.shower2.azimuth
        zeni = self.shower2.zenith
        dire = np.array([np.cos(azi)*np.sin(zeni),
                             np.sin(azi)*np.sin(zeni),
                             np.cos(zeni)])

        # get direction from shower
        direction = self.shower2.direction()

        self.assertEqual(dire.all(), direction.all()) #shall be equal       
        
        print("Finish _direction test")           
    

    
    
    
class Sim_ShowerTest(unittest.TestCase): 
    
    shower3 = SimulatedShower()
    
    def test_sim_shower(self):
        print("Start _simulatedshower test")
        
        # add simulation to shower 
        self.shower3.simulation = "coreas"
        # get simulation
        sim = self.shower3.simulation
        
        self.assertEqual("coreas", sim)
        
        # add Xmax to shower 
        self.shower3.Xmax = 600. *u.g/(u.cm**2)
        # get Xmax
        Xmax = self.shower3.Xmax
        
        self.assertEqual(600. *u.g/(u.cm**2), Xmax)        
        
        
        print("Finish _simulatedshower test")
    


class Reco_ShowerTest(unittest.TestCase): 

    shower4 = ReconstructedShower()
    
        
    def test_recoenergy(self):
        print("Start _recoenergy test")
        recoener = []

        recoener.append(1e18* u.eV )            
        # add energy to shower ---  can only get one primary element per shower
        self.shower4.recoenergy = recoener[-1] 
        # get energy from shower
        recoenergy = self.shower4.recoenergy
        self.assertEqual(recoener[0], recoenergy) #shall be equal       
        
        ## test whether variable can be overwritten              
        #ener.append(3.12e17)            
        ## add energy to shower ---  can only get one primary element per shower
        #self.shower1.add_energy(ener[-1])    
        ## get energy from shower
        #energy = self.shower1.get_energy()        
        #self.assertNotEqual(ener[1], energy) # shall be not equal
                
        print("Finish _recoenergy test")
        
        
    def test_recozenith(self):
        print("Start _recozenith test")
        recozen = []

        recozen.append(70.* u.deg)            
        # add zenithto shower ---  can only get one primary element per shower
        self.shower4.recozenith = recozen[-1]
        # get zenith from shower
        recozenith = self.shower4.recozenith
        self.assertEqual(recozen[0], recozenith) #shall be equal       

        ## test whether variable can be overwritten         
        #zen.append(83.)            
        ## add energy to shower ---  can only get one primary element per shower
        #self.shower1.add_zenith(zen[-1])    
        ## get energy from shower
        #zenith = self.shower1.get_zenith()        
        #self.assertNotEqual(zen[1], zenith) # shall be not equal
                
        print("Finish _recozenith test")
                

    def test_recoazimuth(self):
        print("Start _azimuth test")
        recoazim = []

        recoazim.append(100.* u.deg)            
        # add azimuthto shower ---  can only get one primary element per shower
        self.shower4.recoazimuth = recoazim[-1]    
        # get azimuth from shower
        recoazimuth = self.shower4.recoazimuth
        self.assertEqual(recoazim[0], recoazimuth) #shall be equal       
        
        ## test whether variable can be overwritten         
        #azim.append(185.)            
        ## add energy to shower ---  can only get one primary element per shower
        #self.shower1.add_azimuth(azim[-1])    
        ## get energy from shower
        #azimuth = self.shower1.get_azimuth()        
        #self.assertNotEqual(azim[1], azimuth) # shall be not equal
                
        print("Finish _recoazimuth test")        
                
        

    def test_recoXmax(self):
        print("Start _recoXmax test")
# add Xmax to shower 
        self.shower4.recoXmax = 600. *u.g/(u.cm**2)
        # get Xmax
        recoXmax = self.shower4.recoXmax
        
        self.assertEqual(600. *u.g/(u.cm**2), recoXmax)     
                
        print("Finish _recoXmax test")    
        
        
    shower5 = ReconstructedShower(recoazimuth = 60.* u.deg, recozenith = 45.* u.deg)

    def test_recodirection(self):
        print("Start _recodirection test")
        #shower2 = Shower(azimuth = 60.* u.deg, zenith = 45.* u.deg)
        
        azi = self.shower5.recoazimuth
        zeni = self.shower5.recozenith
        dire = np.array([np.cos(azi)*np.sin(zeni),
                             np.sin(azi)*np.sin(zeni),
                             np.cos(zeni)])

        # get direction from shower
        direction = self.shower5.recodirection()

        self.assertEqual(dire.all(), direction.all()) #shall be equal       
        
        print("Finish _recodirection test")   
       
       
        
        
if __name__ == "__main__":
    unittest.main()
