# Information on scripts
   
   
   
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
