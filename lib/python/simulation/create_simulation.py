
import sys
from sys import argv
import numpy as np
import os
from os.path import split, join, realpath

import matplotlib #.pyplot as plt
#matplotlib.use('Agg')
import pylab as pl
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D


##########################################################################################################
##########################################################################################################

def GRANDtoZHAireS(zen_DANTON=None, azim_DANTON=0):
     """ Convert coordinates from DANTON convention to ZHAireS convention """

     zen = 180. - zen_DANTON
     azim = azim_DANTON - 180.
     if azim>360:
        azim = azim-360.
     elif azim<0.:
        azim = azim+360.
     return [zen,azim]

##########################################################################################################

def array_display(ANTENNAS=None,datamap=None,title=None, core=None):
    import matplotlib.pyplot as pl
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from mpl_toolkits.mplot3d import Axes3D 
    from matplotlib.colors import ListedColormap
    
    if len(ANTENNAS[:,0])!=0:
        fig1 = pl.figure(1,figsize=(5*3.13,3.8*3.13))
        binmap = ListedColormap(['white', 'black'], 'indexed')
        dar=(np.max(ANTENNAS[:,0])-np.min(ANTENNAS[:,0]))/(np.max(ANTENNAS[:,1])-np.min(ANTENNAS[:,1]))
    if dar==0:
          dar=1
    xlbl='X [m]'
    ylbl='Y [m]'
    zlbl='Z [m]'

    ax = pl.gca(projection='3d')
    ax.scatter(ANTENNAS[:,0]*1.,ANTENNAS[:,1],ANTENNAS[:,2],c=datamap)
    ax.scatter(core[0],core[1],core[2],c='red')
    ax.set_title(title)
    ax.view_init(25,-130)
    pl.xlabel(xlbl)
    pl.ylabel(ylbl)
    ax.set_zlabel(zlbl)
    ax.view_init(elev=90., azim=-90.)
    #pl.gca().set_aspect(1,adjustable='box')
    #pl.gca().set_aspect('equal')
    #pl.axis('equal')
    pl.savefig('array.png')

    pl.show()
    return

##########################################################################################################
def reduce_antenna_array(h=None,theta=None,phi=None,ANTENNAS=None,core=[0.,0.,0.],DISPLAY=True):
     """ Reduce the size of the initialized radio array to the shower geometrical footprint by computing the angle between shower and decay-point-to-antenna axes """
     """ theta = zenith in GRAND convention [in deg], h = Xmax height [in m] """
     '''
     Returns:
     --------
     ANTENNAS2: numpy array
        seleted antenna positions nd optional the slopes
     sel: numpy array
        indizes of selected antennas
     '''
     
     
     from os.path import split, join, realpath
     root_dir = realpath(join(split(__file__)[0], "../radio-simus")) # = $PROJECT
     sys.path.append(join(root_dir, "lib", "python"))
     print(join(split(__file__)[0], "../radio-simus"),join(root_dir, "lib", "python") )
     import radio_simus
     from radio_simus.utils import getCerenkovAngle
        
     print('Core = ', core)
     zen,azim = GRANDtoZHAireS(theta,phi)
     zenr = np.radians(zen)
     azimr = np.radians(azim)
     ANTENNAS1 = np.copy(ANTENNAS)

     # TODO fix this with decay position to handle as well neutrinos induced at (0,0,height)

     # Shift antenna array with the randomized core position, core has to be at (0,0) for simulations
     ANTENNAS1[:,0] = ANTENNAS1[:,0]-core[0]
     ANTENNAS1[:,1] = ANTENNAS1[:,1]-core[1]
     ANTENNAS1[:,2] = ANTENNAS1[:,2]#+core[2]

     # Compute angle between shower and decay-point-to-antenna axes
     #u_ant = ANTENNAS1-h*np.array([-np.tan(zenr),0.,1.],dtype=float)
     #u_sh = [np.sin(zenr),0.,-np.cos(zenr)]
     u_ant = ANTENNAS1-np.array([-(h-core[2])*np.tan(zenr)*np.cos(azimr),-(h-core[2])*np.tan(zenr)*np.sin(azimr),h],dtype=float)
     u_ant = (u_ant.T/np.linalg.norm(u_ant,axis=1)).T
     u_sh = np.array([np.cos(azimr)*np.sin(zenr), np.sin(azimr)*np.sin(zenr), -np.cos(zenr)])
     ant_angle = np.rad2deg(np.arccos(np.matmul(u_ant, u_sh)))
     
     # Remove antennas of the initial array that are located outside the "footprint"
     omegar = getCerenkovAngle(h) #[in rad] # Accounting for a footprint four times larger than the Cherenkov angle
     print('--- Cherenkov angle in deg:', omegar, ' times 1.5 for selection' )
     angle_test = ant_angle<=(omegar*1.5)
     sel = np.where(angle_test)[0]
     ANTENNAS2 = ANTENNAS1[sel,:]

     # Remove the farthest antennas to reduce the number of antenna positions to simulate so that this number falls below 1000
     #while np.shape(sel)[0]>9999:
    ##x_ant_max = np.max(np.abs(ANTENNAS2[:,0]))
    ##antisel = np.where(np.abs(ANTENNAS2[:,0])==x_ant_max)[0]
    #r_ant_max = np.max(np.sqrt(ANTENNAS2[:,0]**2+ANTENNAS2[:,1]**2))
    #antisel = np.where(np.sqrt(ANTENNAS2[:,0]**2+ANTENNAS2[:,1]**2)==r_ant_max)[0]
    #ANTENNAS2= np.delete(ANTENNAS2,antisel,0)
    ##sel= np.delete(sel,antisel,0)
    
     ## Shift antenna array with the randomized core position, core has to be at (0,0) for simulations
     #ANTENNAS2[:,0] = ANTENNAS2[:,0]-core[0]
     #ANTENNAS2[:,1] = ANTENNAS2[:,1]-core[1]
     #ANTENNAS2[:,2] = ANTENNAS2[:,2]#+core[2]

     # 3D Display of the radio array
     if DISPLAY:         
        ant_map_i = np.zeros(np.shape(ANTENNAS)[0])
        ant_map_i[sel]=1.
        cc = np.zeros((np.size(ant_map_i),3))
        cc[np.where(ant_map_i==0),:]=[1,1,1]
        #array_display(ANTENNAS,ant_angle,'Shower axis to decay point-antenna axis angle map', core)
        array_display(ANTENNAS,cc,'Selected antenna map', core)

     return ANTENNAS2, sel


##########################################################################################################
def _correct_core(core=[0.,0.,0.], theta=180., azim=0., groundaltitude=0.):
     '''
     Corrects the core from Zhaires inputfile to be corrected for a core height at ground altitude
     
     Parameters:
     core: list
        core coordinates from zhaires inp file
     theta: float
        theta in GRAND conv, deg
     azim: float
        azimuth in GRAND conv, deg
     groundaltitude: float
        observerlevel wrt sealevel
    
     Returns:
     core: list
        core coordinates for hit a groundaltitude
     '''
     
     theta = np.deg2rad(180. - theta) #to zhaires convention
     azim = np.deg2rad(180. + azim)


     h= (core[2]-groundaltitude) # height difference to correct for
     a =np.array([np.sin(theta)*np.cos(azim), np.sin(theta)*np.sin(azim), np.cos(theta)]) # shower direction
     a=a/np.linalg.norm(a)
     a_op = np.array([a[0], a[1],0]) #ortogonal projection on ground to shift for shower axis
     a_op=a_op/np.linalg.norm(a_op)
     d= h* np.tan(theta)# shift length to respect shower axis

     return np.array([core[0],core[1],groundaltitude]) + d*a_op
##########################################################################################################
##########################################################################################################

def create_input_coreas(datadir, fileno, shower, simulation, pos=[]):#, REDUCTION = True):
    ''' creates CoREAS input files
        
        Arguments:
        ----------
        datadir: str
            path to event folder
        fileno: str
            ID number of event
        shower:  dict
            contains infos on shower parameters
        antenna file 
            optional, file containing basic array layout
        simulation: dict
            contains info for simulation
        pos: numpy array
            posx, posy, posz in m, alpha, beta in deg, ID of antenna
            
           
        Returns:
        --------
        0
        
        Creates:
        --------
        Inp, reas, list files to run CoREAS simulations and a general info file for later analysis
        
    '''
    from os.path import split, join, realpath
    root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
    sys.path.append(join(root_dir, "lib", "python"))
    from radio_simus.__init__ import site, Bcoreas, obs_height
    if obs_height != float(shower["core"][2]):
        print("ATTENTION: Check core position vs observer level")

   
    
    primary=None 
    if shower["primary"] == 'proton':
        primary = '14'
    if shower["primary"] == 'iron':
        primary = '5626'
         
    ####### CREATE INP file
    inp_name=datadir+ '/inp/SIM'+fileno+'.inp'
    file= open(inp_name, 'w')

    # work in CoREAS/CORSIKA format/coordinates. direction of propagation: az=0 propagating North
    seed1 = 3*int(fileno)
    seed2 = seed1+1
    seed3 = seed1+2
    
    file.write('RUNNR    {0}\n'.format(fileno))
    file.write('EVTNR    1\n')
    file.write('NSHOW    1\n')# parallel
    file.write('SEED     {0}    0    0\n'.format(seed1)) 
    file.write('SEED     {0}    0    0\n'.format(seed2)) 
    file.write('SEED     {0}    0    0\n'.format(seed3)) 
    file.write('ERANGE    {0:8.2e}    {1:8.2e}\n'.format(shower["energy"]*1e-9, shower["energy"]*1e-9)) # eV in GeV     
    file.write('PRMPAR    {0}\n'.format(primary))
    file.write('THETAP    {0:6.2f}    {1:6.2f}\n'.format(180-shower["zenith"], 180-shower["zenith"])) #GRAND to CoREAS
    file.write('PHIP    {0:6.2f}    {1:6.2f}\n'.format(shower["azimuth"], shower["azimuth"])) 
    #file.write('ECUTS    3.00e-01     3.00e-01     4.010e-04    4.010e-04\n')
    file.write('ECUTS    3.00e-01     1.0     4.010e-04    4.010e-04\n')
    file.write('ELMFLG  T    T\n')
    file.write('THIN     {0:4.2e}     {1:4.2e}     0.00e+00\n'.format(simulation["thinlevel"], shower["energy"]*1e-9 *simulation["thinlevel"]))
    file.write('THINH    1.00e+00     1.00e+02\n')
    file.write('OBSLEV    {0}e2\n'.format(shower["core"][2])) # =groundaltitude
    file.write('ECTMAP  1.00e+11\n')
    file.write('STEPFC  1.00\n')
    file.write('MUMULT  T\n')
    file.write('MUADDI  T\n')
    file.write('HILOW    100.\n')
    file.write('MAXPRT  1\n')
    file.write('MAGNET    {0}    {1}    ({2})\n'.format(Bcoreas[1], Bcoreas[2], site))
    file.write('LONGI    T    5.  T    T\n')
    file.write('RADNKG  5.e5\n')
    file.write('OUTFILE {0}/{1}/{2}/DAT{3}.firstint\n'.format(simulation["path"], simulation["run"], fileno,fileno))
    file.write('DIRECT  {0}/{1}/{2}/\n'.format(simulation["path"], simulation["run"], fileno))
    file.write('EXIT\n')
    
    file.close()


    ####### CREATE REAS file
    reas_name=datadir+'SIM'+fileno+'.reas'
    file2= open(reas_name, 'w')
    
    file2.write('# parameters setting up the spatial observer configuration:\n')
    file2.write('CoreCoordinateNorth = 0      ; in cm\n') # we leave it to (0,0), but instead move the observer positions
    file2.write('CoreCoordinateWest = 0        ; in cm\n')
    file2.write('CoreCoordinateVertical = {0}e2      ; in cm\n'.format(shower["core"][2])) # =groundaltitude
    file2.write('\n')
    file2.write('# parameters setting up the temporal observer configuration:\n')
    file2.write('TimeLowerBoundary = -1 ; in s, only if AutomaticTimeBoundaries set to 0\n')
    file2.write('TimeUpperBoundary = 1 ; in s, only if AutomaticTimeBoundaries set to 0\n')
    file2.write('TimeResolution = {0} ; in s\n'.format(simulation["timebinning"]*1e-9)) #ns to s
    file2.write('ResolutionReductionScale = 0 ; avoid that timing changes\n')
    file2.write('AutomaticTimeBoundaries = 4.e-07 ; 0: off, x: automatic boundaries with width x in s\n')
    file2.write('GroundLevelRefractiveIndex = 1.000292  ; specify refractive index at 0 m asl\n')
    file2.close()
    

    ####### CREATE LIST file - height wrt sealevel
    # read in GP300 antenna files
    if len(pos)>0:
                
        list_name=datadir+ 'SIM'+fileno+'.list'
        file4= open(list_name, 'w')    
        for i in np.arange(0,len(pos.T[0])):
            name="a"+str(pos[i,5]) # name of antenna file
            file4.write('AntennaPosition =  {0:.1f}  {1:.1f}  {2:.1f}  {3}  gamma  2.0  1.e40\n'.format(100*float(pos.T[0,i]), 100*float(pos.T[1,i]), 100*float(pos.T[2,i]), name)) #m in cm
        file4.close()

        
        ####### CREATE INFO file
        info_name=datadir+'SIM'+fileno+'.info'
        file3= open(info_name, 'w')
        
        file3.write('TASK  {0}\n'.format(shower["task"]))
        file3.write('CORE  {0}  {1}  {2}\n'.format(shower["core"][0], shower["core"][1], shower["core"][2])) #in m
        file3.write('EXEC  {0}\n'.format(simulation["executable"]))
        #try:
            #file3.write('SLOPE  {0}\n'.format(simulation["slope"]))
        #except IOError:
            #print('no slope given')
        
        ## TODO get alpha and beta from Turtle according antenna position
        for i in np.arange(0,len(pos.T[0])):
            alpha = pos[i,3]
            beta = pos[i,4]
            ID="{0}".format(pos[i,5])
            file3.write('ANTENNA  {0}  {1:.1f}  {2:.1f}  {3:.1f}  {4}  {5} \n'.format(ID, float(pos.T[0,i]), float(pos.T[1,i]), float(pos.T[2,i]), float(alpha), float(beta))) #m 

        

        # NOTE in case less than 5 antennas in cone delete simulation inputs
        if len(pos.T[0]) < 5:
            import shutil
            print("\n --- ATTENTION: less than 5 antennas in simulation -> delete folder \n")
            #os.system("rm -f {0}".format(datadir))
            shutil.rmtree(datadir)
        else:
            print(str(len(pos.T[0])) + " will be simulated")
        
        file3.close()
        
    return 0





######################################################################

def create_input_zhaires(datadir, fileno, shower, simulation, pos=[]):
        ''' creates ZHAireS input files
        
        Arguments:
        ----------
        datadir: str
            path to event folder
        fileno: str
            ID number of event
        shower:  dict
            contains infos on shower parameters
        antenna file 
            optional, file containing basic array layout
        simulation: dict
            contains info for simulation
        pos: numpy array
            posx, posy, posz in m, alpha, beta in deg, ID of antenna
            
            
        Returns:
        --------
        0
        
        Creates:
        --------
        Inp file to run ZHAireS simulations and a general info file for later analysis
        
        '''
    
    
    
        from os.path import split, join, realpath
        root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
        sys.path.append(join(root_dir, "lib", "python"))
        from radio_simus.__init__ import site, thetageo, phigeo, longitude, latitude, obs_height, Bzhaires
    
    
        ####### CREATE INP file
        inp_name=datadir+ '/SIM'+fileno+'.inp'
        file= open(inp_name, 'w')
        print("save input as: ", inp_name)
        
        primary=None 
        if shower["primary"] == 'proton' or shower["primary"] == 'Proton':
            primary = 'Proton'
        if shower["primary"] == 'iron' or shower["primary"] == 'Iron':
            primary = 'Iron'
        if primary==None:
            print("ATTENTION: primary not yet defined")

        import random
        seed=random.uniform(0, 1)

        task='TaskName '+str(fileno)+ '\n'
        file.write(task)
        prim='PrimaryParticle '+str(primary) + '\n'
        file.write(prim)
        file.write('PrimaryEnergy {0} eV\n'.format(shower["energy"]))
        file.write('PrimaryZenAngle {0:2.1f} deg\n'.format(180-shower["zenith"])) # GRAND to Zhaires
        file.write('PrimaryAzimAngle {0:2.1f} deg Magnetic\n'.format(180+shower["azimuth"]))# GRAND to Zhaires
        file.write('RandomSeed {0:1.3f}\n'.format(seed)) # seems to be  not needed
        #file.write('# mountain slope {0:2.1f} deg\n'.format(alpha*radtodeg))
        file.write('\n')
        file.write('PropagatePrimary On\n')
        file.write('\n')
        file.write('GeomagneticField On\n')
        file.write('\n')
        file.write('SetGlobal RASPASSTimeShift 0.00\n')
        file.write('SetGlobal RASPASSDistance 0.00\n')
        file.write('\n')
        file.write('Thinning   {0:4.2e}  Relative #Thinning is adjusted as per ZHAireS recommendations\n'.format(simulation["thinlevel"]))
        file.write('ThinningWFactor 0.06\n')
        file.write('\n')
        file.write('ElectronCutEnergy 10 MeV\n') # should be discussed
        file.write('ElectronRoughCut 10 MeV\n')
        file.write('GammaCutEnergy 10 MeV\n')
        file.write('GammaRoughCut 10 MeV\n')
        file.write('\n')
        file.write('ForceLowEDecay Never\n')
        file.write('ForceLowEAnnihilation Never\n')
        file.write('\n')
        #file.write('#Location is Ulastai (TREND site 42circ55N 86.68E). where Bgeo is the same as the one used in your GRAND sims:\n')
        #file.write('#phigeo = 182.66; (ie pointing 2.66 degrees East from full North)\n')
        #file.write('#thetageo = 27.05; (pointing down)\n')
        #file.write('#bgeo = 56.5; %muT\n')
        #file.write('#The ground altitude is 2650m.\n')
        file.write('#Linsley\n')
        file.write('Atmosphere 1\n')
  
        file.write('AddSite {0} {1:.1f} deg {2:.1f} deg {3:.1f} m\n'.format(site, latitude, longitude, shower["core"][2]))
        file.write('Site {0}\n'.format(site))
        file.write('Date 1985 10 26\n')
        file.write('GeomagneticField On\n')
        file.write('GeomagneticField {0} uT {1} deg {2} deg\n'.format(str(site), float(Bzhaires[1]), float(Bzhaires[2]), float(Bzhaires[3])))
        file.write('GroundAltitude {0} m\n'.format(shower["core"][2])) # observer level and site height does not have to be core[2]

        file.write('\n')
        file.write('# Save PerShowerData\n')
        file.write('# Parameters of output files.\n')
        file.write('# you wont be needing ground or longitudinal particles\n')
        file.write('\n')
        file.write('PerShowerData Full\n')
        file.write('SaveNotInFile lgtpcles All\n')
        file.write('SaveNotInFile grdpcles All\n')
        file.write('RLimsFile grdpcles 0.000001 m 10 km\n')
        file.write('ResamplingRatio 100\n')
        file.write('RLimsTables 10 m 10 km\n')
        file.write('ELimsTables 10 MeV 1 TeV\n')
        file.write('ObservingLevels 500  0 m 25000 m\n')
        file.write('\n')
        file.write('# All other parameters will be assigned a default value if not set.\n')
        file.write('TotalShowers 1\n')
        file.write('RunsPerProcess Infinite\n')
        file.write('ShowersPerRun 1\n')
        file.write('\n')
        file.write('# In RASPASS all tables are along the shower axis, using planes perpendicular to the axis\n')
        file.write('# The table X axis is in km, or in g/cm2 if Opt a is given.\n')
        file.write('# Opt s supresses header, for gnuplot paste command to work propperly.\n')
        file.write('# Current limitation is that it only works in run time. If tables are exported using AiresExport,\n')
        file.write('# the table content will be ok but the X axis will be wrong\n')
        file.write('\n')
        file.write('#Tables\n')
        file.write('All particles\n')
        file.write('ExportTable 1293 Opt s\n')
        file.write('ExportTable 1293 Opt as\n')
        file.write('ExportTable 1793 Opt s\n')
        file.write('ExportTable 1793 Opt as\n')
        file.write('ExportTable 7293 Opt s\n')
        file.write('ExportTable 7793 Opt s\n')
        file.write('ExportTable 7993 Opt s\n')
        file.write('\n')
        file.write('#e+/e-\n')
        file.write('ExportTable 1205 Opt s\n')
        file.write('ExportTable 1205 Opt as\n')
        file.write('ExportTable 1705 Opt s\n')
        file.write('ExportTable 1705 Opt as\n')
        file.write('ExportTable 7205 Opt s\n')
        file.write('ExportTable 7705 Opt s\n')
        file.write('ExportTable 7905 Opt s\n')
        file.write('\n')
        file.write('#gammas\n')
        file.write('ExportTable 1001 Opt s\n')
        file.write('ExportTable 1001 Opt as\n')
        file.write('ExportTable 1501 Opt s\n')
        file.write('ExportTable 1501 Opt as\n')
        file.write('ExportTable 7001 Opt s\n')
        file.write('ExportTable 7501 Opt s\n')
        file.write('ExportTable 7801 Opt s\n')
        file.write('\n')
        file.write('#mu+/mu-\n')
        file.write('ExportTable 1207 Opt s\n')
        file.write('ExportTable 1207 Opt as\n')
        file.write('ExportTable 1707 Opt s\n')
        file.write('ExportTable 1707 Opt as\n')
        file.write('ExportTable 7207 Opt s\n')
        file.write('ExportTable 7707 Opt s\n')
        file.write('ExportTable 7907 Opt s\n')
        file.write('\n')
        file.write('#pi+/pi-\n')
        file.write('ExportTable 1211 Opt s\n')
        file.write('ExportTable 1211 Opt as\n')
        file.write('ExportTable 1711 Opt s\n')
        file.write('ExportTable 1711 Opt as\n')
        file.write('\n')
        file.write('#k+/k-\n')
        file.write('ExportTable 1213 Opt s\n')
        file.write('ExportTable 1213 Opt as\n')
        file.write('ExportTable 1713 Opt s\n')
        file.write('ExportTable 1713 Opt as\n')
        file.write('\n')
        file.write('#neutrons\n')
        file.write('ExportTable 1021 Opt s\n')
        file.write('ExportTable 1021 Opt as\n')
        file.write('ExportTable 1521 Opt s\n')
        file.write('ExportTable 1521 Opt as\n')
        file.write('\n')
        file.write('#protons\n')
        file.write('ExportTable 1022 Opt s\n')
        file.write('ExportTable 1022 Opt as\n')
        file.write('ExportTable 1522 Opt s\n')
        file.write('ExportTable 1522 Opt as\n')
        file.write('\n')
        file.write('#antiprotons\n')
        file.write('ExportTable 1023 Opt s\n')
        file.write('ExportTable 1023 Opt as\n')
        file.write('ExportTable 1523 Opt s\n')
        file.write('ExportTable 1523 Opt as\n')
        file.write('\n')
        file.write('#nuclei\n')
        file.write('ExportTable 1041 Opt s\n')
        file.write('ExportTable 1041 Opt as\n')
        file.write('ExportTable 1541 Opt s\n')
        file.write('ExportTable 1541 Opt as\n')
        file.write('\n')
        file.write('#och\n')
        file.write('ExportTable 1591 Opt s\n')
        file.write('ExportTable 1591 Opt as\n')
        file.write('ExportTable 7091 Opt s\n')
        file.write('ExportTable 7591 Opt s\n')
        file.write('ExportTable 7891 Opt s\n')
        file.write('\n')
        file.write('#on\n')
        file.write('ExportTable 1592 Opt s\n')
        file.write('ExportTable 1592 Opt as\n')
        file.write('ExportTable 7092 Opt s\n')
        file.write('ExportTable 7592 Opt s\n')
        file.write('ExportTable 7892 Opt s\n')
        file.write('#\n')
        file.write('# ZHAireS v0.28r22\n')
        file.write('#\n')
        file.write('\n')
        file.write('#ZHAiresS Input\n')
        file.write('#ZHAireS On/Off (Default:Off)\n')
        file.write('#Enable the calculation of the electric field. as per the ZHS formalism\n')
        file.write('ZHAireS On\n')
        file.write('\n')
        file.write('#FresnelFreq On/Off (Default:Off)\n')
        file.write('#Perform the calculation in the frequency domain. using the fresnel aproximation\n')
        file.write('FresnelFreq Off\n')
        file.write('\n')
        file.write('#FresnelTime On (Default:Off)\n')
        file.write('#Perform the calculation in the time domain. using the fresnel aproximation\n')
        file.write('FresnelTime On\n')
        file.write('\n')
        file.write('#RefractionIndex n (Default:1.000325. Variable) (note that we are not using the default on purpose)\n')
        file.write('# Refraction Index in EVA?\n')
        file.write('#RefractionIndex 1.00035\n')
        file.write('#ConstRefrIndex\n')
        file.write('\n')
        file.write('#TimeDomainBin t (Default:0.5 ns)\n')
        file.write('#The widht of the bin to be used in time domain calculations\n')
        file.write('TimeDomainBin {0} ns\n'.format(simulation["timebinning"]))
        file.write('AntennaTimeMin -100 ns\n')
        file.write('AntennaTimeMax 850 ns\n')
        file.write('\n')
        file.write('\n')
        file.write('#\n')
        file.write('# End of the fixed part\n')
        file.write('# \n')
        file.write('# This input file has been generated using ProduceInputFile of the AJM suite.\n')
        file.write('# \n')
        
        if len(pos) > 0:
            for i in np.arange(len(pos.T[0])):
                file.write("AddAntenna  {0:11.3f} {1:11.3f} {2:11.3f}\n".format(float(pos[i,0]),float(pos[i,1]),float(pos[i,2])))

        file.close()
        
        
        
        ####### CREATE INFO file
        if len(pos)>0:
                    
            ####### CREATE INFO file
            info_name=datadir+'SIM'+fileno+'.info'
            file3= open(info_name, 'w')
            
            file3.write('TASK  {0}\n'.format(shower["task"]))
            file3.write('CORE  {0}  {1}  {2}\n'.format(shower["core"][0], shower["core"][1], shower["core"][2])) #in m
            file3.write('EXEC  {0}\n'.format(simulation["executable"]))
            #try:
                #file3.write('SLOPE  {0}\n'.format(simulation["slope"]))
            #except IOError:
                #print('no slope given')
            
            ## TODO get alpha and beta from Turtle according antenna position
            for i in np.arange(0,len(pos.T[0])):
                alpha = pos[i,3]
                beta = pos[i,4]
                ID="{0}".format(pos[i,5])
                file3.write('ANTENNA  {0}  {1:.1f}  {2:.1f}  {3:.1f}  {4}  {5} \n'.format(ID, float(pos.T[0,i]), float(pos.T[1,i]), float(pos.T[2,i]), float(alpha), float(beta))) #m 

            

            # NOTE in case less than 5 antennas in cone delete simulation inputs
            if len(pos.T[0]) < 5:
                import shutil
                print("\n --- ATTENTION: less than 5 antennas in simulation -> delete folder \n")
                #os.system("rm -f {0}".format(datadir))
                shutil.rmtree(datadir)
            else:
                print(str(len(pos.T[0])) + " will be simulated")
                
            file3.close()
        else: 
            import shutil
            print("\n --- ATTENTION: no antenna seletected for simulation -> delete folder \n")
            #os.system("rm -f {0}".format(datadir))
            shutil.rmtree(datadir)
            
        
        return 0
   
      
      
##########################################################################################################

### also in frame.py --  to be removed
def get_rotation(zen, az, phigeo=0.72, thetageo=147.43): # included in frame
    """Utility function for getting the rotation matrix between frames
    
    Arguments:
    ----------
    zen: float
        zenith of shower in deg (GRAND)
    az: float
        azimuth of shower in deg (GRAND)
    phigeo: float
        angle magnetic field in deg
    thetageo: float
        angle magnetic field in deg
        
    Returns:
    --------
    numpy array
        rotation matrix
    """
    #magnetic field vector
    s = np.sin(thetageo)
    B = np.array([np.cos(phigeo) * s, np.sin(phigeo) * s,
                     np.cos(thetageo)])
    
    # shower vector   
    s = np.sin(zen)
    v = np.array([np.cos(az) * s, np.sin(az) * s, np.cos(zen)])


    vxB = np.cross(v, B)
    vxB /= np.linalg.norm(vxB)
    vxvxB = np.cross(v, vxB)
    vxvxB /= np.linalg.norm(vxvxB)
    
    return np.array((v, vxB, vxvxB))

##########################################################################################################

def _create_starshape(zen, az, phigeo=0.72, thetageo=147.43, gdalt = 2734, stepsize = 25 ): # Bfield lenghu
    ''' 
    calculation done in CORSIKA coordinates
    '''
    
    zen= 180.-zen # GRAND to CORSIKA
    zen_rad = np.deg2rad(zen)
    az_rad = np.deg2rad(az)
    rot = get_rotation(zen_rad, az_rad, phigeo=phigeo, thetageo=thetageo)
    v = rot[0]
    vxB = rot[1]
    vxvxB = rot[2]
    
    # 160 Antennen
    ang_split= 8.
    #stepsize = 25. #m
    rings= 21

    pos=[]
    for i in np.arange(1,rings):
      for j in np.arange(ang_split):
          xyz = i*(stepsize)*(np.cos(j*(2./(ang_split))*np.pi)*vxB+np.sin(j *(2./(ang_split))*np.pi)*vxvxB)
          #c = xyz[2]/v[2]
          #pos.append([-(xyz[0]-c*v[0]), (xyz[1]-c*v[1]),  gdalt] )
          pos.append([(xyz[0]), xyz[1],  xyz[2]] )

    
    return np.array(pos) 

#########################################################################################################

def _project_starshape(azimuth, zenith, dist_fromxmax, n, core=np.array([0.,0.,0.]) ,max_ang=2.5, thetageo= 147.43, phigeo=0.72):
    
    
    #define shower vector, normal vector plane , core 
    az_rad=np.deg2rad(180.+azimuth)#Note ZHAIRES units used
    zen_rad=np.deg2rad(180.-zenith)


    # shower vector
    v = np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])
    v = v/np.linalg.norm(v)
    
    ### setting starshape
    max_ang = np.deg2rad(2.5) # Most distant antenans are 2degs from axis 
    d1 = dist_fromxmax*np.tan(max_ang)
    step = d1/20. 
    
    xyz0 = _create_starshape(zenith, azimuth, phigeo=0.72, thetageo=147.43, gdalt = 2734, stepsize = step )
    number= len(xyz0.T[0])
    #### star shape pattern in xyz, projected on plane
    xyz=np.zeros([number,3]) 
    rings=int(number/8)
    for i in np.arange(1,rings+1):  
        for j in np.arange(8):
            # projection of shower mountain plane
            b=-np.dot(n,xyz0[(i-1)*8+j])/ np.dot(n, v)
            xyz[(i-1)*8+j]=xyz0[(i-1)*8+j] +b*v + core  #rojected
            
    return xyz
       

          
