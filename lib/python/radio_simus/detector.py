'''
    TODO: several missing functions to be implemented, marked in the code
'''


import numpy as np 
from .__init__ import site, latitude, longitude  ## not yet needed

print("\n================\n Detector set up at:", site,"\n================ \n")

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
        
                
        
        
    ## missing
    # def get_slope() # from position, or antID
    
    # get_position() # get the whole detector
    
    #set_simulated(self, antID)
        ##empty array [antID, 0]
        ##if antID in list of simulated antenna position set 0 to one, return list of antIDs with 1s
        
    #def set_selection(parameters):
        ##empty array [antID, 0]
        ## under any reqirement fulfilled set 0 to one, return list of antIDs with 1s
        
    # def set_antennatype():
    # string
    
    # def get_antennatype():
    # string
    
    # def set_antennaheight():
    # not necessarily the same for all antennas, dependend on antenna type
    
    # def get_antennaheight():
    # not necessarily the same for all antennas, dependend on antenna type
    
#=================================================        

def create_from_file(self, array_file): 
    ''' reading in whole antenna array as file: antID, position and slope
    Arguments:
    ----------
    self
    array_file: str
        filepath containing ID, positions and slope
    
    Returns:
    -------
    -        
    '''
    print("Creating array from: ", array_file, "\n")
    ant_array = np.loadtxt(array_file,  comments="#")
    for i in range(0,len(ant_array.T[0])):
        self.add_position(int(ant_array[i,0]), ant_array[i, 1:4])
        self.add_slope(int(ant_array[i,0]), ant_array[i, 4:])

def get_array(self):
    '''returns array of antenna ID and position for detector, equivalent to get_positions(det) in doc-file
    Arguments:
    ----------
    self
    
    Returns:
    --------
    numpy array
        array of IDs [0] and positions [1:]
    '''
    return np.array(self.array)
    
def get_slopes(self):
    '''returns array of antenna ID and slopes for detector
    Arguments:
    ----------
    self
    
    Returns:
    --------
    numpy array
        array of IDs [0] and slopes [1:]: alpha and beta
    '''
    return np.array(self.slopes)
    

def find_antennaposition(self, ID):
    ''' find antenna position per ID in detector array
    Arguments:
    ----------
    det: object
    ID: int
        desired antenna ID
        
    Returns:
    --------
    positions: array
        desired antenna position
    '''
    array = get_array(self)
    positions = array[:,1:] 
    antIDs = array.T[0]
    index = np.where(antIDs == int(ID))[0][0]
    #print(index)
    return positions[index]


def find_antennaslope(self, ID):
    ''' find antenna slope per ID in detector array
    Arguments:
    ----------
    det: object
    ID: int
        desired antenna ID
        
    Returns:
    --------
    positions: array
        desired antenna slope
    '''
    array = get_slopes(self)
    slopes = array[:,1:] 
    antIDs = array.T[0]
    index = np.where(antIDs == int(ID))[0][0]
    #print(index)
    return slopes[index]



## missing
#def find_antennaID(positions, antIDs, pos):
#    ...
#   return antIDs[index]
