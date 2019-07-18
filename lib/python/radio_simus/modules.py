import os
import sys
import numpy as np
import math
import matplotlib.pyplot as plt



#__all__ = ["compute_ZL", "TopoToAntenna", "_getXmax", "_dist_decay_Xmax"]


#============================================================================


def _geomagnetic_angle(zen,az):
    '''
    Parameters
    ----------
    zen : float 
        zenith of the shower in GRAND, in deg
    az : float 
        azimuth of the shower in GRAND, in deg

    Returns
    -------
    geo_angle: float
        geomagnetic angle in deg
    '''
  
    # Direction of shower
    az_rad = np.degtorad(az)
    zen_rad= np.degtorad(zen)
    v = np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])

    # Direction where Earth mag field points to at Ulastai
    # North + Magnetic North, pointing towards East in this case 
    phigeo = 2.72
    thetageo = 152.95
    print("ATTENTION: theta and phi geo hardcoded: (",thetageo, phigeo ,") deg"
    az_B = np.degtorad(phigeo)
    # Direction where we point to . Inclination=63.05deg
    zen_B = np.degtorad(thetageo)
    B = np.array([np.cos(az_B)*np.sin(zen_B),np.sin(az_B)*np.sin(zen_B),np.cos(zen_B)])

    geo_angle = np.arccos(np.matmul(v, B))
    return np.degtorad(geo_angle)

#============================================================================

def _getXmax(primarytype, energy):
    ''' Xmax value in g/cm2
    Arguments:
    ----------
    primarytype: str
        nature of primary: electron, pion, proton, iron
    energy: float
        primary energy in eV
    
    Returns:
    ----------
    Xmax: float
        Xmax value in g/cm2
    
    Note: factor approximated from https://pos.sissa.it/301/1100/pdf
    '''
       
    print(primarytype, energy)   
    #    current version for tau decays    
    if primarytype=='electron': # aprroximated by gamma shower
        a=82.5 # g/cm2
        c=342.5 #g/cm2
    if primarytype=='pion':  # pion, kaon .... aprroximated by proton
        a=62.5 # g/cm2
        c=357.5 #g/cm2

    if primarytype=='proton'or primarytype=='Proton': # pion, kaon .... approximated by proton
        a=62.5 # g/cm2
        c=357.5 #g/cm2 
        print("ATTENTION: zenith to be corrected for CR")
    if primarytype=='iron' or primarytype=='Iron': # aprroximated by proton
        a=60 # g/cm2 # just approximated 
        c=177.5 #g/cm2
        print("ATTENTION: zenith to be corrected for CR")
    
    Xmax= a*np.log10(energy*10**-12.)+c # eV* 10**-12. to be in TeV

    return Xmax#/abs(np.cos(np.pi-zen2)) # TODO: how to correct for slanted shower
#============================================================================

def _dist_decay_Xmax(zen, injh2, Xmax_primary): 
    ''' Calculate the height of Xmax and the distance injection point to Xmax along the shower axis
    
    Arguments:
    ----------
    zen: float
        GRAND zenith in deg
    injh2: float
        injectionheight above sealevel in m
    Xmax_primary: float
        Xmax in g/cm2 
        
    Returns:
    --------
    h: float
        vertical Xmax_height in m
    ai: float
        Xmax_distance injection to Xmax along shower axis in m
    '''
    
    #% Using isothermal Model
    rho_0 = 1.225*0.001#; % kg/m3 to 0.001g/cm3: 1g/cm3=1000kg/m3, since X given in g/cm2
    M = 0.028966#;  %kg/mol - 1000g/mol
    g = 9.81#; %ms-2
    T = 288.#; % K
    R = 8.32#; J/K/mol , J=kg m2/s2

    zen2 = np.deg2rad(zen)
    
    hD=injh2
    step=10 #m
    if hD>10000:
        step=100 #m
    Xmax_primary= Xmax_primary#* 10. # g/cm2 to kg/m2: 1g/cm2 = 10kg/m2
    gamma=np.pi-zen2 # counterpart of where it goes to
    Re= 6370949 # m, Earth radius
    X=0.
    i=0.
    h=hD
    ai=0
    while X< Xmax_primary:
        i=i+1
        ai=i*step #100. #m
        hi= -Re+np.sqrt(Re**2. + ai**2. + hD**2. + 2.*Re*hD - 2*ai*np.cos(gamma) *(Re+hD))## cos(gamma)= + to - at 90dg
        deltah= abs(h-hi) #(h_i-1 - hi)= delta h
        h=hi # new height
        X=X+ rho_0*np.exp(-g*M*hi/(R*T)) * step*100. #(deltah*100) *abs(1./np.cos(np.pi-zen2)) # Xmax in g/cm2, slanted = Xmax, vertical/ cos(theta); density in g/cm3, h: m->100cm, np.pi-zen2 since it is defined as where the showers comes from, abs(cosine) so correct for minus values
        
    return h, ai # Xmax_height in m, Xmax_distance in m
#============================================================================

def _get_XmaxPosition(primary, energy, zen, azim, injh):
    ''' Calculates vector to Xmax position
    
    Arguments:
    ----------
    primary: str
        primary type, electron or pion
    energy: float
        primary energy in eV
    zen: float
        GRAND zenith in deg
    azim: float
        GRAND azimuth in deg
    injh: float
        injection height wrt sealevel in m
        
    Returns:
    ---------
    new: numpy array
        position of Xmax in m
        
    Note: Only working for neutrinos so far
    '''
    print("ATTENTION: works currently only for neutrinos")
    Xmax_primary = _getXmax(primary, energy)
    # height of xmax, distance decay to xmax
    h, ai = _dist_decay_Xmax(zen, injh, Xmax_primary)
    zenr = np.deg2rad(zen)
    azimr = np.deg2rad(azim)
    u_sh = np.array([np.cos(azimr)*np.sin(zenr), np.sin(azimr)*np.sin(zenr), np.cos(zenr)])
    if primary=="electron" or primary == "pion":
       new = float(ai)*u_sh + np.array([0.,0.,injh ])
    return new

#============================================================================

def _get_CRzenith(zen,injh,GdAlt):
    ''' Corrects the zenith angle for CR respecting Earth curvature, zenith seen by observer
        ---fix for CR (zenith computed @ shower core position
    
    Arguments:
    ----------
    zen: float
        GRAND zenith in deg
    injh: float
        injection height wrt to sealevel in m
    GdAlt: float
        ground altitude of array/observer in m (should be substituted)
    
    Returns:
    --------
    zen_inj: float
        GRAND zenith computed at shower core position in deg
        
    Note: To be included in other functions   
    '''

    #Note: To be included in other functions
    Re= 6370949 # m, Earth radius

    a = np.sqrt((Re + injh)**2. - (Re+GdAlt)**2 *np.sin(np.pi-np.deg2rad(zen))**2) - (Re+GdAlt)*np.cos(np.pi-np.deg2rad(zen))
    zen_inj = np.rad2deg(np.pi-np.arccos((a**2 +(Re+injh)**2 -Re**2)/(2*a*(Re+injh))))
    
    return zen_inj
