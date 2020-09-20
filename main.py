from datetime import datetime
from package.utility import TypeLatLong, TypeOE
from package.orbital import Simulation, OrbitCalculate


"""
TISAT 1                 
1 36799U 10035E   20215.46271709  .00000102  00000-0  16380-4 0  9994
2 36799  98.1977  27.5125 0012384 172.8977 187.2418 14.91391657 545784
"""

if __name__ == "__main__":
    OE = TypeOE(eccentricity=0.0012384, inclination=98.1977, right_ascension=27.5125, argument_of_perigee=172.8977, mean_anomaly=187.2418, mean_motion=14.91391657)
