###
### NOTE: AZ corrected read-in of traces, before sometimes trace.T needed to feed-in
###       added a new didgitization function (other one didnt work for some reason
###

import numpy as np
from scipy.signal import butter, lfilter, resample
from scipy.fftpack import rfft, irfft, rfftfreq


#__all__ = ["add_noise", "digitization", "filter", "run"]

##########################################################################
##########################################################################


def add_noise(vrms, voltages):
    """Add normal random noise on voltages
    inputs : (voltage noise rms, voltages)
    outputs : noisy voltages (time in ns)
    """
    noisy_voltages = np.copy(voltages)
    noisy_voltages[:,1:] = voltages[:,1:] + \
        np.random.normal(0, vrms, size=np.shape(voltages[:,1:]))
    return noisy_voltages

#===========================================================================================================

def digitization(voltages, tsampling):
    """Digitize the voltages at an specific sampling
    inputs : (voltages, sampling rate)
    outputs : digitized voltages
    """
    dt = round(np.mean(np.diff(voltages.T[0]))) #voltages.T[0, 1] - voltages.T[0, 0]
    num = len(voltages.T[0, :]) * int(tsampling / dt)
    if tsampling % dt != 0:
        raise ValueError("Sampling time not multiple of simulation time")
    t_sampled = resample(voltages.T[0, :], num)
    Vx_sampled = resample(voltages.T[1, :], num)
    Vy_sampled = resample(voltages.T[2, :], num)
    Vz_sampled = resample(voltages.T[3, :], num)
    return np.array([t_sampled, Vx_sampled, Vy_sampled, Vz_sampled]).T

#===========================================================================================================

def Digitization_2(v,TSAMPLING):
    """Digitize the voltages at an specific sampling -- v2
    inputs : (voltages, sampling rate)
    outputs : digitized voltages (time in ns)
    """
    v=v.T
    tstep = 1/round(1/(v[0,1] - v[0,0])) # tweak the sh**
    ratio=int(round(TSAMPLING/tstep))
    SAMPLESIZE = int(len(v[0])/ratio)
    vx=np.zeros(SAMPLESIZE)
    vy=np.zeros(SAMPLESIZE)
    vz=np.zeros(SAMPLESIZE)
    tf=np.zeros(SAMPLESIZE)  
    ind=np.arange(0,SAMPLESIZE)*ratio

    if len(ind)>SAMPLESIZE:
        ind=ind[0:TSAMPLING]
    vx[0:len(ind)]=v[1,ind]
    vy[0:len(ind)]=v[2,ind]
    vz[0:len(ind)]=v[3,ind]
    tf[0:len(ind)]=v[0,ind]
    for k in range(len(ind),SAMPLESIZE):
        tf[k]=tf[k-1]+TSAMPLING
    return np.array([tf, vx, vy, vz]).T

#===========================================================================================================

# Filter the voltages at a bandwidth
################################################################
def _butter_bandpass_filter(data, lowcut, highcut, fs):
    """subfunction of filt
    """
    b, a = butter(5, [lowcut / (0.5 * fs), highcut / (0.5 * fs)],
                  btype='band')  # (order, [low, high], btype)
    return lfilter(b, a, data)

#===========================================================================================================

def filters(voltages, FREQMIN=50.e6, FREQMAX=200.e6):
  """ Filter signal v(t) in given bandwidth 
  Parameters
  ----------
   : voltages
      The array of time (s) + voltage (muV) vectors to be filtered
   : FREQMIN 
      The minimal frequency of the bandpass filter (Hz)
   : FREQMAX: 
      The maximal frequency of the bandpass filter (Hz)
      
  Returns
  -------
    numpy array
        time in ns, Voltages (x,y,z)
  Raises
  ------
  Notes
  -----
  At present Butterworth filter only is implemented
  Examples
  --------
  ```
  >>> from signal_treatment import _butter_bandpass_filter
  ```
  """
  
  print("ATTENTION: output traces inversed now")
  t = voltages.T[0]
  # check whether time in s or ns and correct for it
  if t[1]-t[0] > 0.1:
      t*=1e-9 # ns to s
  v = np.array(voltages.T[1:, :])  # Raw signal

  #fs = 1 / np.mean(np.diff(t))  # Compute frequency step
  fs = round(1 / np.mean(np.diff(t)))  # Compute frequency step
  print("Trace sampling frequency: ",fs/1e6,"MHz")
  nCh = np.shape(v.T)[1]
  vout = np.zeros(shape=(len(t), nCh))
  res = t
  for i in range(nCh):
        vi = v[i,:]
        #vout[:, i] = _butter_bandpass_filter(vi, FREQMIN, FREQMAX, fs)
        res = np.append(res,_butter_bandpass_filter(vi, FREQMIN, FREQMAX, fs))
  
  res = np.reshape(res,(nCh+1,len(t)))  # Put it back inright format
  res[0]*=1e9 # s to ns
  return res.T

#===========================================================================================================

##########################################################################
#### RUN THE FULL CHAIN
##########################################################################
### HARDCODED - TODO to be substituted
tsampling = 2 # ns
Vrms = 28 #muV before filtering - NOTE: should be substituted by function returning the value

#===========================================================================================================
def run(efield, zenith_sim, azimuth_sim, alpha_sim=0., beta_sim=0., DISPLAY=1):
        ''' 
        Do the full chain once:
        1. READ IN THE SIMULATED ELECTRIC FIELD TRACE (at higher level) at hand over as parameter
        2. APPLY ANTENNA RESPONSE
        3. ADD STATIONARY NOISE (GALACTIC AND GROUND), VRMS(50-200MHz)= 15muV
        4. FILTER THE TRACE TO THE 50-200MHz WINDOW
        5. DIGITIZATION -- 2ns 
        
        TODO: make it modular so that people can pick the steps they need
    
        
        Arguments:
        ----------
        efield: np array
            electric field trace
        zenith_sim: float
            zenith of shower in deg (GRAND)
        azimuth_sim: float
            azimuth of shower in deg (GRAND)
        alpha, beta: float
            antenna angles, optional, in deg
        DISPLAY: 0,1
            Plotting option off/on
        
        Returns:
        ---------
        trace:
            voltage trace, numpy array: time in ns, voltages (x,y,z)
        
        '''
        import matplotlib.pyplot as plt
        
        print("TSampling in ns: ", tsampling, " , Vrms in muV: ", Vrms )

        ### 2. APPLY ANTENNA RESPONSE
        from radio_simus.computevoltage import get_voltage, compute_antennaresponse
        trace = compute_antennaresponse(efield, zenith_sim, azimuth_sim, alpha=alpha_sim, beta=beta_sim )
        
        #### 2b. deconvolve antenna response - still ongoing work
        #from invert_computevoltage import compute_electicfield
        #electric = compute_electicfield(trace, zenith_sim, azimuth_sim, alpha=alpha_sim, beta=beta_sim )
        
        if DISPLAY==1:
            	    
            plt.figure(1,  facecolor='w', edgecolor='k')
            plt.subplot(311)
            plt.plot(efield.T[0],efield.T[2], label="Ey = EW")
            plt.plot(efield.T[0],efield.T[1], label="Ex = NS")
            plt.plot(efield.T[0],efield.T[3], label="Ez = UP")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Electric field (muV/m)')
            plt.legend(loc='best')
            plt.subplot(312)
            plt.plot(trace.T[0],trace.T[2], label="EW")
            plt.plot(trace.T[0],trace.T[1], label="NS")
            plt.plot(trace.T[0],trace.T[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage (muV)')
            plt.legend(loc='best')
            #plt.subplot(313)
            #plt.plot(electric[0],electric[2], label="Ey = EW")
            #plt.plot(electric[0],electric[1], label="Ex = NS")
            #plt.plot(electric[0],electric[3], label="Ez = UP")
            #plt.xlabel('Time (nsec)')
            #plt.ylabel('Electric field (muV/m)')
            #plt.legend(loc='best')
            
            plt.show()
            #plt.savefig('voltage_antennaresponse.png', bbox_inches='tight')
            

        
        ### 3. ADD STATIONARY NOISE (GALACTIC AND GROUND), VRMS(50-200MHz)= 15muV
        #Vrms = 28 #muV before filtering - NOTE: should be substituted by function returning the value
        trace = add_noise(Vrms, trace) # remove tranposed in signal_treatment
        
        if DISPLAY==1:
            	    
            plt.figure(2,  facecolor='w', edgecolor='k')
            plt.plot(trace.T[0],trace.T[1], label="EW")
            plt.plot(trace.T[0],trace.T[2], label="NS")
            plt.plot(trace.T[0],trace.T[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage + Noise (muV)')
            plt.legend(loc='best')
            
            plt.show()
            #plt.savefig('voltage_addnoise.png', bbox_inches='tight')

            
            
        ### 4. FILTER THE TRACE TO THE 50-200MHz WINDOW
        trace = filters(trace, FREQMIN=50.e6, FREQMAX=200.e6)
    
        
        if DISPLAY==1:
            	    
            plt.figure(3,  facecolor='w', edgecolor='k')
            plt.subplot(211) # time domain
            plt.plot(trace.T[0],trace.T[1], label="EW")
            plt.plot(trace.T[0],trace.T[2], label="NS")
            plt.plot(trace.T[0],trace.T[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage + Noise (muV) - filtered')
            plt.legend(loc='best')
            
            plt.subplot(212) # frequency domain            
            freqs  = np.fft.rfftfreq(len(trace.T[1]),trace.T[0,1]-trace.T[0,0])
            plt.plot(freqs,np.abs(np.fft.rfft(trace.T[1])), label="EW")
            plt.plot(freqs,np.abs(np.fft.rfft(trace.T[2])), label="NS")
            plt.plot(freqs,np.abs(np.fft.rfft(trace.T[3])), label="Vertical")
            plt.xlabel('frequency [Hz]')
            plt.ylabel('Amplitude')
            plt.legend(loc='best')

            plt.show()
            #plt.savefig('voltage_filters.png', bbox_inches='tight')


        ### 5. DIGITIZATION -- 2ns 
        ##trace = digitization(trace,tsampling)
        #trace = Digitization_2(trace,tsampling)
    
        if DISPLAY==1:
            	    
            plt.figure(4,  facecolor='w', edgecolor='k')
            plt.plot(trace.T[0],trace.T[1], label="EW")
            plt.plot(trace.T[0],trace.T[2], label="NS")
            plt.plot(trace.T[0],trace.T[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage + Noise (muV) - digitized')
            plt.legend(loc='best')
            
            plt.show()
            #plt.savefig('voltage_digitisation.png', bbox_inches='tight')
            

        return trace
            
            
#===========================================================================================================
#===========================================================================================================
