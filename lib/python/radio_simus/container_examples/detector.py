from grand import *
import astropy.units as u

origin = grand.ECEF(x = 0 * u.m, y = 0 * u.m, z = 0 * u.m)

print(origin)


"""
see SLACK with Valentin 26Sept
* Define the array (= the positions and ID of all antennas) once
* Internal conversion from any coordinate szstem to ECEF, see example in grand/tools/geomagnet.py
* input type of antenna coordinates could be a Union[grand.ECEF, grand.LTP], ie. local or global coordinates (--> grand pacakage)
* is the underlying data immutable
* detector[0]= 0th antenna data as namedtuple, Antenna(position=[], type='',...)
* use Final decorator (typing_extension
AntennaArry could be a of a list


ToDo
* add positions of a whole array, in ECEF and m
* ...

"""
