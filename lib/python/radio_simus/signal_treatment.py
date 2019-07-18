import numpy as np


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
