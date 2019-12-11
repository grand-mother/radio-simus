''' WORK IN PROGRESS

    a simple copy of CreateInp_GP300.py
    
    run: python3 CreateInp_GP300_1simu.py

    ATTENTION: this scripts run with an input of zhaires simulations (Matias ones)
    hardcoded: 
    * antenna file
    * groundaltitude
    * antenna positions and observation height treated differently than matias does (here I am using one fixed height)
    * thinninglevel and timebinning
    * simus, executable and simpath
    

    ToDo:
    * Documentation
    * effectful reduction of antenna
    * implemnt matias shadowing and antenna selection
    * read in info from config file
    * get the up-to-date antenna file
    *...
'''




import sys
from sys import argv
import numpy as np
import os
from os.path import split, join, realpath


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.mplot3d import Axes3D 


import create_simulation
from create_simulation import _project_starshape, reduce_antenna_array, create_input_coreas, create_input_zhaires

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))

#from radio_simus.modules import _getXmax, _dist_decay_Xmax, _get_CRzenith
#from radio_simus.in_out import inputfromtxt
#from radio_simus.utils import getCerenkovAngle
from radio_simus.modules import _getXmax, _get_CRzenith, _dist_decay_Xmax
########################



# antenna array: file and groundaltitude wrt sealevel
antenna_file = './regular_array_slopes.txt'
alpha=0
beta=0
groundaltitude=2734#.054 #m, from GP300 Zhaires simus

# define name of that specific run - will create a folder with that name
run='test'
# path where to store files
wkdir="./"+run


# TODO: read from general config-file
# Simulation specifics: path for storage on cluster, simulation program, executable, thinninlevel, timebinning
#simus = "zhaires"
#executable = 'zhaires_go'
#simpath = 'cca_at_lyon'


simus="coreas"
executable = 'corsika76400Linux_SIBYLL_gheisha_thin_coreas'
simpath = '/home/fh1-project-huepra/xb6511/coreas/simulations'

# fixed values
thinlevel=1.e-4
timebinning=1. #ns

# reduces the antenna list to the one with expected signals
REDUCTION = True # Reduce antenna array from file


# create only one simulation and handover hardcoded input
ONE=True
showerID="000001"
primarytype="proton"
energy = 1e18 #eV
zen = 180. - 75. # deg (GRAND frame)
az = 180. + 0. # deg (GRAND frame)
injh = 1e6 #m -- not needed for CR showers
task =  "test1" # identification wrt t eg Danton output
core = np.array([0.,0., groundaltitude]) # m , core position wrt to array center


if ONE:
        # set up folders for input file storage
        datadir=wkdir+'/'+showerID+'/'
        if not os.path.exists(datadir):
            #os.makedirs(datadir)
            os.makedirs(datadir+'/inp/')
    
        ########################
        # load shower info from inp file via dictionary
        ########################
        shower = {
          "ID" : showerID,         # shower ID, number of simulation
          "primary" : primarytype,        # primary (electron, pion)
          "energy" : energy,         # eV
          "zenith" : zen,         # deg (GRAND frame)
          "azimuth" : az,          # deg (GRAND frame)
          "injection_height" : injh,    # m (injection height in the local coordinate system)
          "task" : task,    # Identification
          "core" : core    # m, numpy array, core position
          }
        
        simulation = {
            "run" : run, 
            "simulation" : simus, # coreas or zhaires
            "executable" : executable, 
            "path" : simpath,
            "thinlevel" : thinlevel,
            "timebinning" : timebinning
            }
        
        
        
        
        ######## CREATE A GRID ARRAY
        #### get positions of array
        ANTENNAS = np.loadtxt(antenna_file,comments="#") # should contain posx, posy, posz, alpha, beta
       
        
        ### TODO ADD RETRO TO CALCULATE ALPHA BETA IF NOT HANDED OVER
        print('hardcoded alpha and beta in antenna file')

        
        ##### ALL ANTENNAS
        ANT  = []
        
        if REDUCTION:
        #### get parameters to select antennas by cone
            Xmax_primary= _getXmax(shower["primary"], shower["energy"]) # eV
            # corrected zen needed for CR and inclined showers
            if shower["primary"] == "proton" or shower["primary"] == "proton" or shower["primary"] == "Iron" or shower["primary"] == "iron":
                zen2 = _get_CRzenith(shower["zenith"],shower["injection_height"],shower["core"][2])
            else:
                zen2 = np.copy(shower["zenith"])
                #print(zen2)
            Xmax_height, Xmax_distance = _dist_decay_Xmax(zen2, shower["injection_height"], Xmax_primary) 
            #print("Xmax_height, Xmax_distance: ", Xmax_height, Xmax_distance)

            ### Reduce the radio array to the shower geometrical footprint (we account for a footprint twice larger than the Cherenkov angle)
            ### Antenna positions corrected for shower core, als
            pos, sel = reduce_antenna_array(Xmax_height,shower["zenith"],shower["azimuth"],ANTENNAS[:,1:4],shower["core"],DISPLAY=True) # at this point positions in m
            print('---- number of antennas: ', len(sel), ' out of ', len(ANTENNAS.T[0]))
            
            SHOW=False
            if SHOW:
                core_show= np.array([0.,0.,groundaltitude])
                shower_dir=np.zeros([200,3])
                # TODO: substitude by get_direction
                a =np.array([np.sin(np.deg2rad(180-shower["zenith"]))*np.cos(np.deg2rad(shower["azimuth"]+180)), np.sin(np.deg2rad(180-shower["zenith"]))*np.sin(np.deg2rad(shower["azimuth"]+180)), np.cos(np.deg2rad(180-shower["zenith"]))]) # shower direction
                for i in np.arange(0,200):     
                    shower_dir[i]= (i)*100 *a + core
                
                fig2 = plt.figure(2,figsize=(5,3.8))
                ax = plt.gca(projection='3d')
                ax.scatter(ANTENNAS[:,1]*1.,ANTENNAS[:,2],ANTENNAS[:,3])
                ax.scatter(ANTENNAS[sel[:],1]*1.,ANTENNAS[sel[:],2],ANTENNAS[sel[:],3])
                ax.scatter(shower_dir[:,0]*1.,shower_dir[:,1],shower_dir[:,2])
                plt.show()
        else:
            pos = np.copy(ANTENNAS)
            pos[:,0:3] = pos[:,0:3] - np.array([core[0], core[1], 0.]) # correct for shower core
            sel = range(len(pos.T[0])) #TODO list of all indizes  
            
        if len(sel) >0:
            for i in range(len(pos.T[0])):
                name="{}".format(sel[i]) #index of selected antenna 
                ANT.append([pos[i,0], pos[i,1], pos[i,2], alpha, beta, name])

            ANT=np.array(ANT)



        if simus == "coreas":
            create_input_coreas(datadir, showerID,  shower, simulation, ANT)
        if simus == "zhaires":
            create_input_zhaires(datadir, showerID,  shower, simulation, ANT)
          


     

