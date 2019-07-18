

import numpy as np 
#from .__init__ import site, latitude, longitude  ## not yet needed

#print("here", site)

class detector:
    ''' info on antenna array

        
    '''
    def __init__(self):
        self.max_length = 1
        
        self.array = []# attribute references
        self.slopes = []# attribute references


    
    #def add_antID(self, antID): # check whether ID already exist
        #if len(self.antID)==self.max_length:
            #print("Warning for ant ", str(self.antID),": Can't add antID=", str(antID)," -- already set to: ", self.antID[0])
        #elif antID is None:
            #print("Warning: antID not given")
        #else:
            #self.antID.append(antID)
            
    #def get_antID(self):
        #if len(self.antID)==0:
            #print("Warning: no antID set")
        #else:
            #return self.antID[0] # str
        
     
     
    def add_position(self, antID, position):
        # add antID and position to array
        if antID in self.array:
            print("AntID already exists")
        elif antID is None:
            print("Warning: AntID not given")
        elif position in self.array:
            print("position already exists")
        elif position is None:
            print("Warning: position not given")
        else:
            self.array.append([antID, position[0], position[1], position[2]])
        
    #def get_position(self, antID=None, position=None):
    
    
    def add_slope(self, antID, slope):
        # add antID and slope (alpha, beta) to array
        if antID in self.slopes:
            print("AntID already exists")
        elif antID is None:
            print("Warning: AntID not given")
        #elif position in self.array:
            #print("position already exists")
        #elif position is None:
            #print("Warning: position not given")
        else:
            self.slopes.append([antID, slope[0], slope[1]])
        
                
        
        
    
    #def get_slope() # from position, or antID
    
    # get_position() # get the whole detector
    
    #set_simulated(self, antID)
        ##empty array [antID, 0]
        ##if antID in list of simulated antenna position set 0 to one, return list of antIDs with 1s
        
    #def set_selection():
        ##empty array [antID, 0]
        ## under any reqirement fulfilled set 0 to one, return list of antIDs with 1s

#=================================================        

def create_from_file(self, array_file): 
    # reading in whole antenna array as file
    # antID == index of position
    print(array_file)
    ant_array = np.loadtxt(array_file,  comments="#")
    for i in range(0,len(ant_array.T[0])):     
        self.add_position(i, ant_array[i, 0:3])
        self.add_slope(i, ant_array[i, 3:])

def get_array(self):
        # returns array of antenna ID and position
        return np.array(self.array[:][1:])
    
def get_slopes(self):
        # returns array of antenna ID and position
        return np.array(self.slopes[:][1:])
    

