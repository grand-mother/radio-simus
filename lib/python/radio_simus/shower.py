"""
My first try to write a python class

Shower class: shall contain all infos on the shower
"""



class shower:
    ''' info on shower parameter '''
    def __init__(self, name):
        self.name = name # attribute references
        self.zenith = []    # How to limit to one number per list .....
        self.azimuth = [] #instantiation
        
        self.simulation = []
        
        self.recoenergy = []
        
    def add_zenith(self, zenith):
        self.zenith.append(zenith)

    def get_zenith(self):
        return self.zenith
        
    def add_azimuth(self, azimuth):
        self.azimuth.append(azimuth)
        
    def get_azimuth(self):
        return self.azimuth 
        
    def add_all(self, zenith, azimuth):
        self.add_zenith(zenith)
        self.add_azimuth(azimuth)
        
    #def get_energy(self):
        
    #def get_primary(self):
 
 
class sim_shower(shower):
    ''' info on simulations '''
    def add_simulation(self, simulation):
        self.simulation.append(simulation)
        
    def get_simulation(self):
        print(self.simulation)


class reco_shower(shower):
    ''' info on reconstrcuted values '''
    def add_recoenergy(self, recoenergy):
        self.recoenergy.append(recoenergy)
        
    def get_recoenergy(self):
        print(self.recoenergy)        
        
        
#=================================================        


sh1 = shower("1")
#sh1.get_zenith()

sh1.add_all(93., 187)
print(sh1.get_zenith())

sh1.add_azimuth(180)
print(sh1.get_azimuth())


sh2 = sim_shower("2")
sh2.add_simulation("coreas")
sh2.get_simulation()


