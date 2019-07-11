import numpy as np

#===========================================================================================================
def load_trace(directory, index, suffix=".trace"):
    """Load data from a trace file (ascii file)
    
    Arguments:
    ----------
    directory: str
        path to folder
    index: int
        antenna ID
    suffix: str
        file ending
        
    Returns:
    --------
    numpy array
        field trace
        
    Raises:
    --------
    IOError:
        adopt to ID naming
        
    Note: currently only usable for Zhaires simulations    
    """
    try:
        path = "{:}/a{:}{:}".format(directory, index, suffix)
        with open(path, "r") as f:
            return np.array([map(float, line.split()) for line in f])
    except IOError:
        path = "{0}/a{1:04d}{2}".format(directory, index+1, suffix)
        with open(path, "r") as f:
            return np.array([map(float, line.split()) for line in f])
        
#===========================================================================================================

def getn(h):
    """Get the refractive index
    
    Arguments:
    ----------
    h: float
        height in m wrt to sealevel
        
    Returns:
    --------
    float
        refractive index at height h
        
    Note: Reference: Zhaires (see email M. Tueros 25/11/2016)
    """
    # h in meters
    return 1. + 325E-06 * np.exp(-0.1218E-03 * h)

#===========================================================================================================

def getCerenkovAngle(h):
   """Get the Cerenkov angle
    
    Arguments:
    ----------
    h: float
        height in m wrt to sealevel
        
    Returns:
    --------
    float
        Cherenkov angle at height h in deg
   """
   return np.rad2deg(np.arccos(1. / getn(h)))

#===========================================================================================================
def get_integratedn(injh2, position):
    '''Calculates the integrated n from a specific height (0,0,height) to the observer position
        --- works currently only for neutrinos
        
    Arguments:
    ----------
    injh2: float
        injection height wrt sealevel in m (0.,0., injh2)
    position: numpy array
        observer position wrt sealevel in m 
    
    Returns:
    --------
    n : float
        integrated refractive index along a path from injection to observer
    
    Note: assumption coordinate system so that tau decay at (0.,0, injectionheight)
    Note: calculation of integrated n implemented similar as in Zhaires
    '''
    
    Re= 6370949 # m, Earth radius
    ########
    # line of sight
    ux= position[0] -0.
    uy= position[1] -0.
    uz= position[2] -injh2
    
    nint= 10000 # so many sub tracks alng line of sight
    # vector k along the line of sight
    kx=ux/nint
    ky=uy/nint
    kz=uz/nint
    
    #current positions, start with injh as first point of emission
    currpx=0.
    currpy=0.
    currpz=injh2
    currh=currpz # just in that case here, since particle injected directly induce a shower
    
    ns=325E-06
    kr=-0.1218E-03
    
    #print "inhh, antenna height ", injh2, position[2]
    summe=0.
    for i in range(0,nint):
        nextpx=currpx+kx
        nextpy=currpy+ky
        nextpz=currpz+kz
        
        nextR=np.sqrt( nextpx*nextpx +nextpy*nextpy )
        nexth= ( np.sqrt((( injh2 - nextpz  ) + Re) * (( injh2  - nextpz  ) + Re) + nextR*nextR) - Re) /1e3
        
        if (abs(currh-nexth)>1e-10 ):
            summe=summe+ (  np.exp(kr*nexth) -   np.exp(kr*currh)  )/ (kr*( nexth - currh) )
        else:
            summe=summe+ np.exp(kr*currh)
        
        currpx=nextpx
        currpy=nextpy
        currpz=nextpy
        currR=nextR
        currh=nexth
        
        
    avn= ns*summe/nint
    n= 1.+ avn
    
    return  n # integrated n

#===========================================================================================================
def mag(x):
    ''' Calculates the length of a vector or the absolute value
    '''
    return np.sqrt(x.dot(x))

#===========================================================================================================
def p2p(trace):
    ''' Calculates peak to peak values
    
    Arguments:
    ----------
    trace: np.array with 1 or 4 columns
            signal trace
    Returns:
    --------
    list: floats
        peak-to-peak values
    '''
    if trace.ndim==2: # simply assume np.array([t, x, y, z])
        return max(trace.T[1])-min(trace.T[1]), max(trace.T[2])-min(trace.T[2]), max(trace.T[3])-min(trace.T[3])
    elif trace.ndim==1:
        return max(trace)-min(trace)
    else:
        print("in p2p(): dimensions not correct")

#===========================================================================================================
def hilbert_env(signal):
    ''' 
    Hilbert envelope - abs(analyticcal signal)
    Arguments:
    ----------
    signal: np.array
        signal trace
    Returns:
    --------
    amplitude_envelope: np.array
        Hilbert envelope
        
    '''
    from scipy.signal import hilbert
    analytic_signal = hilbert(signal)
    amplitude_envelope = np.abs(analytic_signal)
    return amplitude_envelope

#===========================================================================================================
def hilbert_peak(time, signal):
    ''' Calculates time and amplitude of peak
    
    Arguments:
    ----------
    time: numpy array
        time in ns
    signal: numpy array
        signal trace in muV or muV/m
        
    Returns:
    --------
    maximum of Hilbert envelope in muV or muV/m: float
    time of maximum in ns : float

    '''
    envelope=hilbert_env(signal)
    #Get time and amp
    
    return max(envelope), time[np.where(signal == signal.max())][0]
 
 
#=========================================================================================================== 
def _getAngle(refpos=[0.,0.,1e6],theta=None,azim=None,ANTENNAS=None, core=[0.,0.,0.]): # theta and azim in Grand convention
    """ Get angle between antenna and shower axis (injection point or Xmax)
        Arguments:
        ----------
        refpos: numpy array
            sofar position of injection or Xmax
        theta: float
            GRAND zenith in deg
        azim: float
            GRAND azimuth in deg
        ANTENNAS: numpy array
            observer position in m
        core: numpy array
            optional, core position in m, not yet used
            
        Returns:
        --------
        float
            Angle in deg between shower axis and vector reference position to observer position
            
        Note: currently only working for neutrinos
        TODO make it working for CRs, currently only for neutrinos
    """

    zenr = np.radians(theta)
    azimr= np.radians(azim)
    ANTENNAS1 = np.copy(ANTENNAS)

    # Compute angle between shower and decay-point-to-antenna axes
    u_ant = ANTENNAS1-refpos
    u_ant = (u_ant/np.linalg.norm(u_ant))

    u_sh = [np.cos(azimr)*np.sin(zenr), np.sin(azimr)*np.sin(zenr), np.cos(zenr)]
    ant_angle = np.arccos(np.matmul(u_ant, u_sh))

    return np.rad2deg(ant_angle)
