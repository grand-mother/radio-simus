import numpy as np

#===========================================================================================================
def load_trace(directory, index, suffix=".trace"):
    """Load data from a trace file
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

       Reference:
        Zhaires (see email M. Tueros 25/11/2016)
    """
    # h in meters
    return 1. + 325E-06 * np.exp(-0.1218E-03 * h)

def getCerenkovAngle(h):
   """Get the Cerenkov angle
   """
   return np.arccos(1. / getn(h))
#===========================================================================================================
def get_integratedn(zen2, injh2, position):
    
    # assumption coordinate system so that tau decay at (0.,0, injectionheight)
    # calculation of integrated n implemented similar as in Zhaires
    
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
    return np.sqrt(x.dot(x))

#===========================================================================================================
def p2p(trace):
    ''' 
    Parameters:
        trace: np.array with 1 or 4 columns
    Returns:
        peak-to-peak values: floats
    '''
    if trace.ndim==2:
        return max(trace.T[1])-min(trace.T[1]), max(trace.T[2])-min(trace.T[2]), max(trace.T[3])-min(trace.T[3])
    elif trace.ndim==1:
        return max(trace)-min(trace)
    else:
        print("in p2p(): dimensions not correct")

#===========================================================================================================
def hilbert_env(signal):
    ''' 
    Hilbert envelope
    Parameters
        signal: np.array
    Returns:
        amplitude_envelope: np.array
        
    '''
    from scipy.signal import hilbert
    analytic_signal = hilbert(signal)
    amplitude_envelope = np.abs(analytic_signal)
    return amplitude_envelope

#===========================================================================================================
def hilbert_peak(signal):
    ''' 
    Returns time and amplitude of peak
    '''
    envelope=hilbert_env(signal)
    #Get time and amp
    return 0
    



