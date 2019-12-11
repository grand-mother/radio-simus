# Information on scripts
Collection of example scripts using the python radio-simus library 
A list of not-fully-ready, but somehow usable scripts:
    
    example_simtohdf5.py: Script loops over event folders and produces hdf5 files per event, containing efield and voltage trace (only tested with coreas simulation so far)
        -- reads in simulations (only coreas tested), stores event info and efield traces in one hdf5 file
        -- performs voltage computation and saves trace in hdf5 file (optional)
        -- performs little analysis on triggering and stores it in hdf5 file (optional)
        -- coomented: some additional comments to use atsropy tables etc
    Usage: python example_simtohdf5.py <path to event folders> <zhaires/coreas>
    
    example_usingclass.py:  Example on how to use classes shower and detector, read-in only hdf5 format, all-antenna files
        -- reads in analysis trigger for events, create a list of events (class objects) and  trigger 1/0 to class attributes
        -- create a png with statistic for triggering
        -- Note: a plot of voltage distribution in in 2d would be nice --> TODO list
    Usage: python3 example_usingclass.py <folder event set>
    
    example_plot_2D.py: Example on how to do a 2D plot, show the radio footprint in the array
        -- read in original array list (from config-file)
        -- plots the p2p distribution componentwise
        -- reads in event hdf5 files
    Usage: python3 example_plot_2D.py <path to event file> <e/v/efield/voltages>
    
    example_plot_trace.py: Script to read in a single antenna trace
        -- Plot trace (PLOT)
        -- Plot Hilbert envelope (HILB)
        -- calculate Angle (ANGLE) -- to be fixed for CR (modules.py)
    Usage: python3 example_plot_trace.py <path to event> <keyword:efield/voltages/...> <antID>
    
    example_voltage_hdf5.py: Example how to read in event info and traces from hdf5 file and run processing of field
        -- reads event info and traces from all-antenna hdf5 file
        -- show how to use the list/options of modules in standard_processing
        -- saves the voltages traces in hdf5 file
    Usage: python3 example_voltage_hdf5.py <path to event folder> 
    
    
   
   
## to example_simtohdf5.py:
Currently, the script can produce one single hdf5-file per antenna, for SINGLE=True. 
The data structure for the produced hdf5-file is:

 **EXAMPLE ALL=True: Structure of hdf5 files, containing all antennas of event [event*.hdf5]:**
 
    /event: stores ID, position and slopes of antennas in event
    /event _meta: example:
        datatype:, 
        - {name: ant_ID, datatype: string}, 
        - {name: pos_x, unit: m, datatype: float64}, 
        - {name: pos_y, unit: m, datatype: float64}, 
        - {name: pos_z, unit: m, datatype: float64}, 
        - {name: alpha, unit: deg, datatype: float64}, 
        - {name: beta, unit: deg, datatype: float64}, 
        meta:, 
        ID: '000002', 
        azimuth: 6.78, 
        core: [15397.397388851985, -530.2640221736006, 2734.0], 
        energy: 4.05e+17, 
        injection_height: 10000000.0, 
        primary: proton, 
        simulation: coreas, 
        task: 36.E.4e17_Z.97_A.7_La.39_Lo.92_H.3078_D.0, 
        zenith: 96.65, 
        
        
    /event/<ant_ID>
    
    /event/<ant_ID>/efield : stores Time, Ex, Ey, Ez as astropy table
    /event/<ant_ID>/efield _meta: example:
        datatype:, 
        - {name: Time, unit: ns, datatype: float64}, 
        - {name: Ex, unit: u V / m, datatype: float64}, 
        - {name: Ey, unit: u V / m, datatype: float64}, 
        - {name: Ez, unit: u V / m, datatype: float64}, 
        meta:, 
        position: [999.9973888519853, -2499.964022173601, 2734.0], 
        slopes: [0.0, 0.0], 
        
    /event/<ant_ID>/voltages: stores Time, Vx, Vy, Vz as astropy table
    /event/<ant_ID>/voltage _meta: example:
        datatype:, 
        - {name: Time, unit: ns, datatype: float64}, 
        - {name: Vx, unit: u V, datatype: float64}, 
        - {name: Vy, unit: u V, datatype: float64}, 
        - {name: Vz, unit: u V, datatype: float64}, 
        meta:, 
        position: [999.9973888519853, -2499.964022173601, 2734.0], 
        slopes: [0.0, 0.0], 
        voltage: antennaresponse, 
        
        
    /analysis:
            trigger: [threshold(aggr) in muV, 1/0 in any, 1/0 in xy, threshold(cons) in muV, 1/0 in any, 1/0 in xy] # 1 = triggered
        p2p: [Ex, Ey, Ez, Exy, Exz] #p2p-value in muV
        
    
## Access to data ( ID == antenna ID ):
    from astropy.table import Table

    Time           Ex                   Ey                    Ez         
    ns         u V / m              u V / m               u V / m       
    ------- -------------------- -------------------- ---------------------
    10668.0                  0.0                  0.0                   0.0
    10669.0                  0.0                  0.0                   0.0
    10670.0                  0.0                  0.0                   0.0
    10671.0                  0.0                  0.0                   0.0
    ...
    11147.0 -0.22002369833168597 -0.28481404837811314  -0.28438298034182286
    11148.0 -0.21189438651797549 -0.24039673025105382 -0.023664534043256736
    11149.0 -0.14237827913618956  -0.1427215146215885   0.04836943976456008
    Length = 482 rows

    
    >>> print(f["Ex"].unit)
    u V / m
    
    >>> print(f.info)
    <Table length=482>
    name  dtype    unit 
    ---- ------- -------
    Time float64      ns
    Ex float64 u V / m
    Ey float64 u V / m
    Ez float64 u V / m

    >>> print(f.meta)
    {'position': [-3500.0, 1000.0, 2734.0], 'slopes': [0.0, 0.0]}
    >>> print(f.meta["position"])
    [-3500.0, 1000.0, 2734.0]

    
    
    
    
## Access event data -- example:
 
    g=Table.read("data.hdf5", path="/event")
    
    >>> print(g.meta)
    {'ID': '000001', 'azimuth': <Quantity 136.38 deg>, 'core': <Quantity [   0.,    0., 2734.] m>, 'energy': <Quantity 2.88e+17 eV>, 'injection_height': <Quantity 10000000. m>, 'primary': 'proton', 'simulation': 'coreas', 'task': '7.E.3e17_Z.96_A.136_La.39_Lo.92_H.3320_D.0', 'zenith': <Quantity 96.26 deg>}
    
    >>> print(g.meta["primary"])
    proton

    
    >>> print(g.info)
    <Table length=128>
    name   dtype  unit
    ------ ------- ----
    ant_ID  bytes3     
    pos_x float64    m
    pos_y float64    m
    pos_z float64    m
    alpha float64  deg
    beta float64  deg

    >>> print(g["alpha"])
    alpha
    deg 
    -----
    0.0
    ...
    0.0
    0.0
    0.0
    Length = 128 rows
    
    >>> print(g["alpha"].unit)
    deg
    
    >>> print(g.meta["energy"])
    2.88e+17 eV
    >>> print(g.meta["energy"].unit)
    eV





        
        
        
<!--        trigger: [threshold(aggr) in muV, 1/0 in any, 1/0 in xy, threshold(cons) in muV, 1/0 in any, 1/0 in xy] # 1 = triggered
        p2p: [Ex, Ey, Ez, Exy, Exz] #p2p-value in muV-->
        
 
        
ToDo: add units in meta info (python dict)) for the event correctly. Currently, the units are set via the read-in from the raw simulation files.
Units:
* ID: '-' (ID of the event)
* azimuth: GRAND deg (azimuth of shower)
* zenith: GRAND deg (zenith of shower)
* core: m (shower core position wrt to array center)
* energy: eV (primary energy)
* injection_height: m (for CR = 100km, for nu = injection height above sealevel) 
* primary: string (primary nature)
* simulation: string (coreas/zhaires/RM,...)
* task: string (task ID, eg from DANTON or RETRO output)
* trigger threshold: muV
* p2p values: muV




## Example of configuration file (still ongoing design phase)
The config files contains all information needed to run simulations as well as for later processing. 

 Given example test.config:

    # define site  ## not yet needed
    SITE  Lenghu
    LONG  92.334037  deg
    LAT  38.870398  deg 
    OBSHEIGHT  2734.0  m

    #antenna file
    ARRAY  /home/laval1NS/zilles/CoREAS/regular_array_slopes.txt

    #magnetic field: 
    THETAGEO  147.43  deg
    PHIGEO  0.72  deg
    B_COREAS  28.17  28.17  # Bx  Bz  uT 
    B_ZHAIRES  54.021  57.43  0.72  ## strengh F in uT, inclination I in deg, declination D in deg

    # definition sigma in muV (50-200MHz)
    VRMS  15
    VRMS2  28  # before filtering, approximated
    TSAMPLING  2  ns  # for digitisation

    #antenna responses
    ANTX  /home/laval1NS/zilles/radio-simus/lib/python/radio_simus/GRAND_antenna/HorizonAntenna_SNarm_leff_loaded.npy
    ANTY  /home/laval1NS/zilles/radio-simus/lib/python/radio_simus/GRAND_antenna/HorizonAntenna_EWarm_leff_loaded.npy
    ANTZ  /home/laval1NS/zilles/radio-simus/lib/python/radio_simus/GRAND_antenna/HorizonAntenna_Zarm_leff_loaded.npy

    # # AddSite Ulastai 42.55 deg 86.68 deg 1500 
    # # GeomagneticField 56.5000 uT 63.18 deg 2.72 deg
    
