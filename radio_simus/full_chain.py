################################
#### by A. Zilles, last update 22 May 2019
################################


#!/usr/bin/env python
import os
from os.path import  join
import sys
#import math
import numpy as np

from scipy.fftpack import rfft, irfft, rfftfreq

import pylab as pl
import matplotlib.pyplot as plt

# from GRAND packages
from signal_treatment import filters, add_noise,digitization, Digitization_2
from computevoltage import get_voltage, compute_antennaresponse


### HARDCODED - TODO to be substituted
tsampling = 2* 1e-9 # ns -> s
Vrms = 28 #muV before filtering - NOTE: should be substituted by function returning the value





#===========================================================================================================
def run(opt_input,path,l, zenith_sim, azimuth_sim, alpha_sim=0., beta_sim=0., DISPLAY=1):
        ''' 
        Do the full chain once:
        1. READ IN THE SIMULATED ELECTRIC FIELD TRACE
        2. APPLY ANTENNA RESPONSE
        3. ADD STATIONARY NOISE (GALACTIC AND GROUND), VRMS(50-200MHz)= 15muV
        4. FILTER THE TRACE TO THE 50-200MHz WINDOW
        5. DIGITIZATION -- 2ns 
        
        TODO: make it modular so that people can pick the steps they need
    
        
        Arguments:
        ----------
        opt_input: str
            manual, txt etc - not yet fixed, needed for read-in
        path: str
            path to folder containing traces
        l: int
            antenna ID
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
            voltage trace, numpy array
        
        '''
        

        ### 1. READ IN THE SIMULATED ELECTRIC FIELD TRACE
        efieldtxt=path+'/a'+str(l)+'.trace'
        print('\n**  Efield file:',efieldtxt)
        
        # ------ LOADING FROM TEXT FILE
        try:
	    # Loadinfg traces
            #time1_sim, Ex_sim,Ey_sim,Ez_sim = np.loadtxt(efieldtxt,usecols=(0,1,2,3),unpack=True)
            efield = np.loadtxt(efieldtxt,usecols=(0,1,2,3),unpack=True)

        except IOError:
            print('IOError: file not found')
        
        # ------- NOTE to be substituted by loading from hdf5 file
        
        
        ### 2. APPLY ANTENNA RESPONSE
        
        trace = compute_antennaresponse(efield, zenith_sim, azimuth_sim, alpha=alpha_sim, beta=beta_sim )
        
        
        #### 2b. deconvolve antenna response - still ongoing work
        #from invert_computevoltage import compute_electicfield
        #electric = compute_electicfield(trace, zenith_sim, azimuth_sim, alpha=alpha_sim, beta=beta_sim )
        

        ####plots
        if DISPLAY==1:
            	    
            plt.figure(1,  facecolor='w', edgecolor='k')
            plt.subplot(311)
            plt.plot(efield[0],efield[2], label="Ey = EW")
            plt.plot(efield[0],efield[1], label="Ex = NS")
            plt.plot(efield[0],efield[3], label="Ez = UP")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Electric field (muV/m)')
            plt.legend(loc='best')
            plt.subplot(312)
            plt.plot(trace[0]*1e9,trace[2], label="EW")
            plt.plot(trace[0]*1e9,trace[1], label="NS")
            plt.plot(trace[0]*1e9,trace[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage (muV)')
            plt.legend(loc='best')
            #plt.subplot(313)
            #plt.plot(electric[0],electric[2], label="Ey = EW")
            #plt.plot(electric[0],electric[1], label="Ex = NS")
            #plt.plot(electric[0],electric[3], label="Ez = UP")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Electric field (muV/m)')
            plt.legend(loc='best')
            
            plt.show()
            #plt.savefig('voltage.png', bbox_inches='tight')
            

        
        
        ### 3. ADD STATIONARY NOISE (GALACTIC AND GROUND), VRMS(50-200MHz)= 15muV
        
        #Vrms = 28 #muV before filtering - NOTE: should be substituted by function returning the value
        trace = add_noise(Vrms, trace) # remove tranposed in signal_treatment
                       

        if DISPLAY==1:
            	    
            plt.figure(2,  facecolor='w', edgecolor='k')
            plt.plot(trace[0]*1e9,trace[1], label="EW")
            plt.plot(trace[0]*1e9,trace[2], label="NS")
            plt.plot(trace[0]*1e9,trace[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage + Noise (muV)')
            plt.legend(loc='best')
            
            plt.show()
            
            
        ### 4. FILTER THE TRACE TO THE 50-200MHz WINDOW
        
        trace = filters(trace, FREQMIN=50.e6, FREQMAX=200.e6)
        
        if DISPLAY==1:
            	    
            plt.figure(3,  facecolor='w', edgecolor='k')
            plt.subplot(211) # time domain
            plt.plot(trace[0]*1e9,trace[1], label="EW")
            plt.plot(trace[0]*1e9,trace[2], label="NS")
            plt.plot(trace[0]*1e9,trace[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage + Noise (muV) - filtered')
            plt.legend(loc='best')
            
            plt.subplot(212) # frequency domain            
            freqs  = np.fft.rfftfreq(trace.shape[1],trace[0,1]-trace[0,0])
            plt.plot(freqs,np.abs(np.fft.rfft(trace[1])), label="EW")
            plt.plot(freqs,np.abs(np.fft.rfft(trace[2])), label="NS")
            plt.plot(freqs,np.abs(np.fft.rfft(trace[3])), label="Vertical")
            plt.xlabel('frequency [Hz]')
            plt.ylabel('Amplitude')
            plt.legend(loc='best')

            plt.show()
        
        ### 5. DIGITIZATION -- 2ns 
        #tsampling = 2* 1e-9 # ns -> s
        #trace = digitization(trace, tsampling) # NOTE does not work for me
        trace =Digitization_2(trace,tsampling)
        
        if DISPLAY==1:
            	    
            plt.figure(4,  facecolor='w', edgecolor='k')
            plt.plot(trace[0]*1e9,trace[1], label="EW")
            plt.plot(trace[0]*1e9,trace[2], label="NS")
            plt.plot(trace[0]*1e9,trace[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage + Noise (muV) - digitized')
            plt.legend(loc='best')
            
            plt.show()

        return trace
            
            
#===========================================================================================================
#===========================================================================================================

if __name__ == '__main__':

    if ( len(sys.argv)<7 ):
        print("""
        This is an example script how to use the simulation chain (run() function) starting with the electric field to the voltage trace detected at the antenna level. 
        -- So far just working for a single antenna.
        
        All angles are to be expressed in degrees and in GRAND convention.
        
        Usage for single antenna:
            python full_chain.py [path to traces] [antID] [zenith] [azimuth] [alpha] [beta]
            
            example: python full_chain.py ./ 58 85 205 10 5 
            
        """)
        sys.exit(0)

    ## Plotting option
    PLOT=1


    #folder containing the trace
    path = sys.argv[1] 
    # Antenna ID
    antID = sys.argv[2]
    
    opt_input = 'manual' # TODO: default for now
    #if opt_input=='txt':
        ## Read the ZHAireS input (.inp) file to extract the primary type, the energy, the injection height and the direction
        #inp_file = str(sys.argv[1])
        #zenith_sim,azimuth_sim = inputfromtxt(inp_file)

    if opt_input=='manual':
        zenith_sim = float(sys.argv[3]) #deg
        azimuth_sim = float(sys.argv[4]) #deg
        alpha = float(sys.argv[5]) #deg
        beta = float(sys.argv[6]) #deg
        
    ### Start the run
    voltage = run(opt_input,path,antID, zenith_sim, azimuth_sim, alpha, beta, PLOT)

    ### save as txt file for later use
    name = "./out_{:}.dat".format(antID)#os.path.join('./', )
    with open(name, "w+") as FILE:
        for i in range(0, len(voltage[0])):
            args = (voltage[0, i], voltage[1, i], voltage[2, i], voltage[3, i])
            try:
                print("%e	%1.3e	%1.3e	%1.3e" % args, end='\n', file=FILE)
            except SyntaxError:
                print >> FILE, "%1.3e	%1.3e	%1.3e	%1.3e" % args
                
                
