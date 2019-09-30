
from .shower import Shower, SimulatedShower

import astropy.units as u


if __name__ == "__main__":


    simulated = SimulatedShower(primary = "tau-", energy = 1E+19 * u.eV,
                                simulation="ZHAires")
    print("This is a simulated shower: ", simulated)
