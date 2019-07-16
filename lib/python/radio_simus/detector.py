class detector:
    ''' info on antenna array

        
    '''
    def __init__(self):
        self.max_length = 1
        
        self.antID = [] # attribute references


    
    def add_antID(self, antID): # check whether ID already exist
        if len(self.antID)==self.max_length:
            print("Warning for ant ", str(self.antID),": Can't add antID=", str(antID)," -- already set to: ", self.antID[0])
        elif antID is None:
            print("Warning: antID not given")
        else:
            self.antID.append(antID)
            
    #def get_antID(self):
        #if len(self.antID)==0:
            #print("Warning: no antID set")
        #else:
            #return self.antID[0] # str
        
        
    #def add_position(antID, position) # assign position to antID
    
    #def add_antenna(antID, position) # add ant_ID and position
    
    #def get_slope() # from position, or antID
    
    # get_position() # get the whole detector
