# Information on scripts
Collection of example scripts using the python radio-simus library 
   
   
   
## to example_simtohdf5.py:
Currently, the script can produce one single hdf5-file per event, for ALL=True. 
The data structure for the produced hdf5-file is:

 **Structure of hdf5 files, containing all antennas of event [event*.hdf5]:**
 
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
    /eventt/<ant_ID>/efield _meta: example:
        datatype:, 
        - {name: Time, unit: ns, datatype: float64}, 
        - {name: Ex, unit: u V / m, datatype: float64}, 
        - {name: Ey, unit: u V / m, datatype: float64}, 
        - {name: Ez, unit: u V / m, datatype: float64}, 
        meta:, 
        position: [999.9973888519853, -2499.964022173601, 2734.0], 
        slopes: [0.0, 0.0], 
        
    /event/<ant_ID>/voltages: stores Time, Vx, Vy, Vz as astropy table
    /eventt/<ant_ID>/efield _meta: example:
        datatype:, 
        - {name: Time, unit: ns, datatype: float64}, 
        - {name: Vx, unit: u V, datatype: float64}, 
        - {name: Vy, unit: u V, datatype: float64}, 
        - {name: Vz, unit: u V, datatype: float64}, 
        meta:, 
        position: [999.9973888519853, -2499.964022173601, 2734.0], 
        slopes: [0.0, 0.0], 
        voltage: antennaresponse, 
        
ToDo: add units in meta info (python dict)) for the event. Currently, the units are set via the read-in from the raw simulation files.
Units:
* ID: - (ID of the event)
* azimuth: GRAND deg (azimuth of shower)
* zenith: GRAND deg (zenith of shower)
* core: m (shower core position wrt to array center)
* energy: eV (primary energy)
* injection_height: m (for CR = 100km, for nu = injection height above sealevel) 
* primary: string (primary nature)
* simulation: string (coreas/zhaires/RM,...)
* task: string (task ID, eg from DANTON or RETRO output)




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
    B_ZHAIRES  54.021  57.43  0.72  #54.021 uT 57.43 deg 0.72 deg

    # definition sigma in muV (50-200MHz)
    VRMS  15
    VRMS2  28  # before filtering
    TSAMPLING  2  ns  # for digitisation

    #antenna responses
    ANTX  /home/laval1NS/zilles/radio-simus/lib/python/radio_simus/GRAND_antenna/HorizonAntenna_SNarm_leff_loaded.npy
    ANTY  /home/laval1NS/zilles/radio-simus/lib/python/radio_simus/GRAND_antenna/HorizonAntenna_EWarm_leff_loaded.npy
    ANTZ  /home/laval1NS/zilles/radio-simus/lib/python/radio_simus/GRAND_antenna/HorizonAntenna_Zarm_leff_loaded.npy

    # # AddSite Ulastai 42.55 deg 86.68 deg 1500 
    # # GeomagneticField 56.5000 uT 63.18 deg 2.72 deg
    
