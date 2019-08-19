# -*- coding: utf-8 -*-
"""
Unit tests for the radio_simus.version module
"""

import unittest
import sys
import numpy as np

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
import radio_simus.shower as sh
#from framework import git


class ShowerTest(unittest.TestCase):
    """Unit tests for the version module"""
    
    #def _init_(self): # Constructor
    shower1 = sh.shower()

    def test_showerID(self):
        print("Start _showerID test")
        ID = []
        for i in range(2):
            shower_ID=None
            ID.append("name"+str(i))            
            # add ID to shower ---  can only get one ID per shower
            self.shower1.add_showerID(str("name"+str(i)))    
            # get ID from shower
            shower_ID = self.shower1.get_showerID()
                        
            if i==0:
                self.assertEqual(ID[0], shower_ID) #shall be equal
                self.assertEqual(len(ID[0]), len(shower_ID))
            if i==1:
                self.assertNotEqual(ID[1], shower_ID) # shall be not equal
                self.assertEqual(len(ID[1]), len(shower_ID))
                
        print("Finish _showerID test")
        

    def test_primary(self):
        print("Start _primary test")
        prim = []

        prim.append("proton")            
        # add primary to shower ---  can only get one primary element per shower
        self.shower1.add_primary(prim[-1])    
        # get primary from shower
        primary = self.shower1.get_primary()
        self.assertEqual(prim[0], primary) #shall be equal
        
               
        prim.append("electron")            
        # add primary to shower ---  can only get one primary element per shower
        self.shower1.add_primary(prim[-1])    
        # get primary from shower
        primary = self.shower1.get_primary()        
        self.assertNotEqual(prim[1], primary) # shall be not equal
                
        print("Finish _primary test")
        
        
    def test_energy(self):
        print("Start _energy test")
        ener = []

        ener.append(1e18)            
        # add energy to shower ---  can only get one primary element per shower
        self.shower1.add_energy(ener[-1])    
        # get energy from shower
        energy = self.shower1.get_energy()
        self.assertEqual(ener[0], energy) #shall be equal       
        
        ener.append(3.12e17)            
        # add energy to shower ---  can only get one primary element per shower
        self.shower1.add_energy(ener[-1])    
        # get energy from shower
        energy = self.shower1.get_energy()        
        self.assertNotEqual(ener[1], energy) # shall be not equal
                
        print("Finish _energy test")
        
        
    def test_zenith(self):
        print("Start _zenith test")
        zen = []

        zen.append(70.)            
        # add zenithto shower ---  can only get one primary element per shower
        self.shower1.add_zenith(zen[-1])    
        # get zenith from shower
        zenith = self.shower1.get_zenith()
        self.assertEqual(zen[0], zenith) #shall be equal       
        
        zen.append(83.)            
        # add energy to shower ---  can only get one primary element per shower
        self.shower1.add_zenith(zen[-1])    
        # get energy from shower
        zenith = self.shower1.get_zenith()        
        self.assertNotEqual(zen[1], zenith) # shall be not equal
                
        print("Finish _zenith test")
                

    def test_azimuth(self):
        print("Start _azimuth test")
        azim = []

        azim.append(100.)            
        # add azimuthto shower ---  can only get one primary element per shower
        self.shower1.add_azimuth(azim[-1])    
        # get azimuth from shower
        azimuth = self.shower1.get_azimuth()
        self.assertEqual(azim[0], azimuth) #shall be equal       
        
        azim.append(185.)            
        # add energy to shower ---  can only get one primary element per shower
        self.shower1.add_azimuth(azim[-1])    
        # get energy from shower
        azimuth = self.shower1.get_azimuth()        
        self.assertNotEqual(azim[1], azimuth) # shall be not equal
                
        print("Finish _azimuth test")        
        
        
    def test_injectionheight(self):
        print("Start _injectionheight test")
        injh = []

        injh.append(1500.)            
        # add injectionheightto shower ---  can only get one primary element per shower
        self.shower1.add_injectionheight(injh[-1])    
        # get injectionheight from shower
        injectionheight = self.shower1.get_injectionheight()
        self.assertEqual(injh[0], injectionheight) #shall be equal       
        
        injh.append(2000.)            
        # add energy to shower ---  can only get one primary element per shower
        self.shower1.add_injectionheight(injh[-1])    
        # get energy from shower
        injectionheight = self.shower1.get_injectionheight()        
        self.assertNotEqual(injh[1], injectionheight) # shall be not equal
                
        print("Finish _injectionheight test")             
        
    
    
    shower2 = sh.shower()
    
    def test_all(self):
        print("Start _all test")
        
        # add all to shower ---  can only get one primary element per shower
        self.shower2.add_all("2", "pion", 3e17, 85., 10, None)    
        # get all from shower
        all_shower = self.shower2.get_all()
        self.assertEqual(3e17, all_shower[2]) #shall be not equal since set before bz other test       
        
        print("Finish _all test")
    
    
    
    
class Sim_ShowerTest(unittest.TestCase): 
    
    shower3 = sh.sim_shower()
    
    def test_sim_shower(self):
        print("Start _sim_shower test")
        
        # add all to shower (inherited) ---  can only get one primary element per shower
        self.shower3.add_all("2", "pion", 3e17, 85., 10, None)  
        # add simulation to shower 
        self.shower3.add_simulation("coreas")
        # get all from shower
        all_shower = self.shower3.get_all()
        # get simulation
        sim = self.shower3.get_simulation()
        
        self.assertEqual(3e17, all_shower[2]) #shall be not equal since set before bz other test       
        self.assertEqual("coreas", sim[0])
        
        print("Finish _sim_shower test")
    


class Reco_ShowerTest(unittest.TestCase): 

    shower4 = sh.reco_shower()
    
    def test_recoenergy(self):
        print("Start _recoenergy test")
        ener = []

        ener.append(1e18)            
        # add energy to shower ---  can only get one primary element per shower
        self.shower4.add_recoenergy(ener[-1])    
        # get energy from shower
        energy = self.shower4.get_recoenergy()
        self.assertEqual(ener[0], energy) #shall be equal       
        
        ener.append(3.12e17)            
        # add energy to shower ---  can only get one primary element per shower
        self.shower4.add_recoenergy(ener[-1])    
        # get energy from shower
        energy = self.shower4.get_recoenergy()        
        self.assertNotEqual(ener[1], energy) # shall be not equal
                
        print("Finish _recoenergy test")
        
        
        
    def test_recoXmax(self):
        print("Start _recoXmax test")
        X = []

        X.append(400.)            
        # add Xmax to shower ---  can only get one primary element per shower
        self.shower4.add_recoXmax(X[-1])    
        # get Xmax from shower
        Xmax = self.shower4.get_recoXmax()
        self.assertEqual(X[0], Xmax) #shall be equal       
        
        X.append(553.)            
        # add Xmax to shower ---  can only get one primary element per shower
        self.shower4.add_recoXmax(X[-1])    
        # get Xmax from shower
        Xmax = self.shower4.get_recoXmax()        
        self.assertNotEqual(X[1], Xmax) # shall be not equal
                
        print("Finish _recoXmax test")     
        
        
    def test_recozenith(self):
        print("Start _recozenith test")
        zen = []

        zen.append(93.)            
        # add zenith to shower ---  can only get one primary element per shower
        self.shower4.add_recozenith(zen[-1])    
        # get zenith from shower
        zenith = self.shower4.get_recozenith()
        self.assertEqual(zen[0], zenith) #shall be equal       
        
        zen.append(85.)            
        # add zenith to shower ---  can only get one primary element per shower
        self.shower4.add_recozenith(zen[-1])    
        # get zenith from shower
        zenith = self.shower4.get_recozenith()        
        self.assertNotEqual(zen[1], zenith) # shall be not equal
                
        print("Finish _recozenith test") 
        

    def test_recoazimuth(self):
        print("Start _recoazimuth test")
        azim = []

        azim.append(153.)            
        # add azimuth to shower ---  can only get one primary element per shower
        self.shower4.add_recoazimuth(azim[-1])    
        # get azimuth from shower
        azimuth = self.shower4.get_recoazimuth()
        self.assertEqual(azim[0], azimuth) #shall be equal       
        
        azim.append(200.)            
        # add azimuth to shower ---  can only get one primary element per shower
        self.shower4.add_recoazimuth(azim[-1])    
        # get azimuth from shower
        azimuth = self.shower4.get_recoazimuth()        
        self.assertNotEqual(azim[1], azimuth) # shall be not equal
                
        print("Finish _recoazimuth test") 
        
        
        
    shower5 = sh.reco_shower()
    
    def test_recoall(self):
        print("Start _recoall test")
        
        # add all to shower ---  can only get one primary element per shower
        self.shower5.add_all("2", "pion", 3e17, 85., 10, 1500)    
        # get all from shower
        all_shower = self.shower5.get_all()
        self.assertEqual(3e17, all_shower[2]) #shall be not equal since set before bz other test       
        
        # add reco to shower ---  can only get one primary element per shower
        self.shower5.add_recoall(655, 4.5e19, 89., 110)    
        # get all from shower
        recoall_shower = self.shower5.get_recoall()
        self.assertEqual(89., recoall_shower[2]) #shall be not equal since set before bz other test
        
        print("Finish _recoall test")
       
       
        
        
if __name__ == "__main__":
    unittest.main()
