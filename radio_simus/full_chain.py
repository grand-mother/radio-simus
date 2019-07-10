################################
#### by A. Zilles
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
from in_out import inputfromtxt


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

        
        ### 2. APPLY ANTENNA RESPONSE
        trace = compute_antennaresponse(efield, zenith_sim, azimuth_sim, alpha=alpha_sim, beta=beta_sim )
        
        #### 2b. deconvolve antenna response - still ongoing work
        #from invert_computevoltage import compute_electicfield
        #electric = compute_electicfield(trace, zenith_sim, azimuth_sim, alpha=alpha_sim, beta=beta_sim )
        
        #print(trace)

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
        #print(trace.T[0])               
        
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
        trace =digitization(trace,tsampling)
        
        if DISPLAY==1:
            	    
            plt.figure(4,  facecolor='w', edgecolor='k')
            plt.plot(trace[0]*1e9,trace[1], label="EW")
            plt.plot(trace[0]*1e9,trace[2], label="NS")
            plt.plot(trace[0]*1e9,trace[3], label="Vertical")
            plt.xlabel('Time (nsec)')
            plt.ylabel('Voltage + Noise (muV) - digitized')
            plt.legend(loc='best')
            
            plt.show()
        print(trace)

        return trace
            
            
#===========================================================================================================
#===========================================================================================================

if __name__ == '__main__':

    if ( len(sys.argv)<1 ):
        print("""
        This is an example script how to use the simulation chain (run() function) starting with the electric field to the voltage trace detected at the antenna level. 
        -- So far just working for a single antenna.
        
        All angles are to be expressed in degrees and in GRAND convention.
        
        Usage for single antenna:
            python full_chain.py [path to traces] [antID] [zenith] [azimuth] [alpha] [beta]
            
            example: python full_chain.py ./ 58 85 205 10 5 (in manual mode)
                     python full_chain.py ./                (in txt mode)
            
        """)
        sys.exit(0)

    ## Plotting option
    PLOT=1

    #folder containing the trace
    path = sys.argv[1] 
    
    ### READ-IN Options, TODO to be set as a parameter
    opt_input = 'manual' 
    
    if opt_input=='txt': # for antenna arrays
        # Read the ZHAireS input (.inp) file to extract the primary type, the energy, the injection height and the direction
        # NOTE alpha and beta are now hardcoded for whole array - to be fixed
        inp_file = str(sys.argv[1])
        showerID=str(inp_file).split("/")[-2] # should be equivalent to folder name
        inputfile = path+'/inp/'+showerID+'.inp'
        zenith_sim,azimuth_sim,energy,injh,primarytype= inputfromtxt(inputfile)
        alpha = 5.
        beta = 0.
        print('Attention: alpha, beta have to be adjusted')
        
        ### define the loop
        # NOTE ZHaires, antenna positions given in m, but here doesnt matter 
        positions=np.genfromtxt(path + '/antpos.dat') # Zhaires
        start=0
        end=len(positions)
        

    if opt_input=='manual': # for single antennas
        antID = int(sys.argv[2])# Antenna ID
        zenith_sim = float(sys.argv[3]) #deg, GRAND conv.
        azimuth_sim = float(sys.argv[4]) #deg, GRAND conv.
        alpha = float(sys.argv[5]) #deg, defining slope
        beta = float(sys.argv[6]) #deg, defining slope
        
        start = antID
        end = antID+1
        

    

    print('Now looping over',end-start,'antenna(s) in folder',path)
    for l in range(start,end):   
        
        ### 1. READ IN THE SIMULATED ELECTRIC FIELD TRACE
            
        # ------- NOTE to be substituted by loading from hdf5 file
        if opt_input=='manual' or 'txt':
            try:
                ## Define the path to the electric field trace
                ## NOTE: Zhaires input: t / ns, Ex / muV/m, Ey / muV/m, Ez / muV/m 
                efieldtxt=path+'/a'+str(l)+'.trace' # ZHaires
                print('\n**  Efield file:',efieldtxt)
                
                # Loading traces
                efield = np.loadtxt(efieldtxt,usecols=(0,1,2,3),unpack=True)
                
                # NOTE: that CoREAS is working in cgs units -> convert electric field in muV/m and time in ns
                #       Antenna positions are given in cm -> convert to m (not needed here) 

            except IOError:
                print('IOError: file not found')       
        
        #if opt_input=='hdf5':# TODO
            #name = path+'/table'+str(l)+'.hdf5'
            #efield = Table.read(name, path=name)
            
            
        ### Start the run
        try:
            voltage = run(efield, zenith_sim, azimuth_sim, alpha, beta, PLOT)


            ### save as txt file for later use in <path>
            # ------- NOTE to be substituted by loading from hdf5 file
            if opt_input=='manual' or 'txt':
                name = path+"/out_{:}.dat".format(l)#os.path.join('./', )
                with open(name, "w+") as FILE:
                    for i in range(0, len(voltage[0])):
                        args = (voltage[0, i], voltage[1, i], voltage[2, i], voltage[3, i])
                        try:
                            FILE.write("{0:e}  {1:f}  {2:f}  {3:f}'\n'".format(voltage[0, i], voltage[1, i], voltage[2, i], voltage[3, i]))
                        except SyntaxError:
                            print >> FILE, "%1.3e	%1.3e	%1.3e	%1.3e" % args
                            
                
            
                    
        except IndexError:
            print('--- ATTENTION: no voltage trace computed or stored')
            continue      
