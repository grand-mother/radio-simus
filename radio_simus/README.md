
# Current status and targeted goals - transfer to a common data format

## Output format of simulations - historical view :)

* ZHAireS/RASPASS: 
    - output file: all traces are stored in one txt file (units: s, V/m). For a simplier use, we split-up the traces into single-antenna files (units: ns, muV/m).
    **Be aware of** that for RASPASS the atmosphere gets "flipped". Therefore, the z-component of the antenna position in meters in the output file as well as the z-component of the electric field have to be 'back-flipped' (multiply by -1)!
    - angle convention :
    theta_zhaires= 180deg - theta_grand
    phi_zhaires= phi_grand + 180 deg
    - fix by Matias: thanks to Matias, the big output file is now directly split-up into separate electric-field files for each antenna 'a<i>.dat' (ns, muV/m) and a separate file listing the antenna positions 'antpos.dat'
    
* CoREAS:
    - output file: the traces for the individual antenna positions are saved in individual textfiles (<name><i>.txt, **units: cgs!** ) with names defined in the antenna-list file beforehand. Multiply the electric field by 2.99792458e4 to transfer the electric field from cgs to SI units (V/m).  The antenna positions in centimeters are handed over as input as a separated file (*.list). 
    - angle convention: the angle are defined by the direction the shower travels to. 
    Azimuth= 0 deg (magnetic North), Azimuth=90 deg (West) -> same as in GRAND
    Zenith = between particle momentum and **negative** z-axis!
    theta_coreas= 180deg -theta_grand
    phi_coreas= phi_grand
    
    theta_coreas= theta_zhaires
    phi_coreas= phi_zhaires - 180deg

* Calculate the voltage traces:
    The module 'compute_voltage' accepts numpy.array with time in s, Ex, Ey,Ez in muV/m. The shower direction has to be defined in GRAND conventions. (Helper function 'compute_antennaresponse' accepts only ns)
    It returns  the voltage traces as a numpy.array with time in s, Vx,Vy,Vz in muV
    
    
## ToDos:
* module reading in CoREAS shower parameters (using astropy.unit)
* module reading in CoREAS antenna list (using astropy.unit)
* module reading in ZHAireS shower parameters (using astropy.unit)
* module reading in ZHAireS antenna list (using astropy.unit)

* module storing traces in hdf5 format (using astropy.unit)
* how to store shower parameters? -> not framework
* how to store antenna informations? -> not framework

* adopt full_chain to hdf5 format
* add astropy logger
* add unit tests


