# Local frame transforms for pulse shape computations.
import numpy
from .__init__ import phigeo, thetageo

def get_rotation(zen, az, phigeo=phigeo, thetageo=thetageo):
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
    s = numpy.sin(thetageo)
    B = numpy.array([numpy.cos(phigeo) * s, numpy.sin(phigeo) * s,
                     numpy.cos(thetageo)])
    
    # shower vector   
    s = numpy.sin(zen)
    v = numpy.array([numpy.cos(az) * s, numpy.sin(az) * s, numpy.cos(zen)])


    vxB = numpy.cross(v, B)
    vxB /= numpy.linalg.norm(vxB)
    vxvxB = numpy.cross(v, vxB)
    vxvxB /= numpy.linalg.norm(vxvxB)
    
    return numpy.array((v, vxB, vxvxB))

# ---------------------------------------------------------

def UVWGetter(cx=0., cy=0., cz=0., zen, az, phigeo=phigeo, thetageo=thetageo):
    """Closure for getting coordinates in the shower frame.
    
    Arguments:
    ----------
    cx, cy, cz: float
        center eg. of antenna positions
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
        rotated vector in shower frame
    """
    R = get_rotation(zen, az, phigeo, thetageo)
    origin = numpy.array((cx, cy, cz))

    def GetUVW(pos):
       return numpy.dot(R, pos - origin)
    return GetUVW

# ---------------------------------------------------------

def XYZGetter(cx, cy, cz, zen, az, phigeo=phigeo, thetageo=thetageo):
    """Closure for getting back to the main frame
        
    Arguments:
    ----------
    cx, cy, cz: float
        center eg. of antenna positions
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
        rotated vector from shower frame in main frame
    """

    Rt = get_rotation(zen, az, phigeo, thetageo).T
    origin = numpy.array((cx, cy, cz))

    def GetXYZ(pos):
        return numpy.dot(Rt, pos) + origin
    return GetXYZ


