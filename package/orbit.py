import numpy as np
from datetime import datetime, timedelta
from typing import Tuple
from utility import (CONSTANT, TypeOE, TypeLatLong, TypeXYZ, normalize, calMeanAnomaly, calEccentricAnomaly, calVecInPQW, toECIfromLatLong, getRie, getMatECItoPQW)
from orbital_elements import calc_oe_from_sv


class Simulation(object):
    def __init__(self, OE: TypeOE, observer: TypeLatLong, E0: float = 1, time: datetime = datetime.now()) -> None:
        self.__OE = OE.copy()
        self.eccentric_anomaly = calEccentricAnomaly(OE['mean_anomaly'], OE['eccentricity'], E0)
        self.__observer = observer  
        self.__time = time  

    def _update(self, dt: float) -> None:
        """Recalculate MeanAbnomaly  and EccentricAnomaly after dt
        """
        self.__OE['mean_anomaly'] = calMeanAnomaly(self.__OE['mean_anomaly'], self.__OE['semimajor_axis'], dt)
        self.eccentric_anomaly = calEccentricAnomaly(self.__OE['mean_anomaly'], self.__OE['eccentricity'], self.eccentric_anomaly)

    def _calPosPQW(self, dt: float) -> np.array:
        """Calculate Position Vector of Satellite in PQW Frame
        """
        self._update(dt)
        return calVecInPQW(self.__OE['semimajor_axis'], self.__OE['eccentricity'], self.eccentric_anomaly)

    def _getObserverECI(self, dt: float) -> np.array:
        return toECIfromLatLong(self.__observer, time=self.__time + timedelta(seconds=dt))

    def _getPosInECI(self, dt:float) -> np.array:
        pos_pqw = self._calPosPQW(dt)
        toECI = getMatECItoPQW(i=self.__OE['inclination'], omega=self.__OE['argument_of_perigee'], sigma=self.__OE['right_ascension'])
        coord = toECI.T @ pos_pqw
        return coord

    def getRie(self, dt: float) -> np.array:
        """Calculate Rie from self.__time + dt(second)
        """
        return getRie(self.__time + timedelta(seconds=dt))

    def getAllCoords(self, dt: float) -> Tuple[np.array, np.array]:
        pos = self._getPosInECI(dt)
        observer = self._getObserverECI(dt)
        return observer, pos


class OrbitCalculate(object):

    def __init__(self, sim: Simulation) -> None:
        self.sim = sim

    def calculate(self, time_points: Tuple[float, float, float], r0: float=600, repeated: int=10000, alias: float = 0.1) -> TypeOE:

        if (repeated < 0 or alias < 0):
            raise ValueError("Repeated must be int. Repeated and Alias must be positive.")
        observers = []
        directions = []

        for dt in time_points:
            observer, pos = self.sim.getAllCoords(dt)
            direction = normalize(pos - observer)
            observers.append(observer)
            directions.append(direction)
        print(directions)

        t1 = time_points[0] - time_points[1]
        t3 = time_points[2] - time_points[1]
        t = t3 - t1

        D0 = directions[0] @ np.cross(directions[1], directions[2])
        D11 = observers[0] @ np.cross(directions[1], directions[2])
        D21 = observers[1] @ np.cross(directions[1], directions[2])
        D31 = observers[2] @ np.cross(directions[1], directions[2])

        D12 = observers[0] @ np.cross(directions[0], directions[2])
        D22 = observers[1] @ np.cross(directions[0], directions[2])
        D32 = observers[2] @ np.cross(directions[0], directions[2])

        D13 = observers[0] @ np.cross(directions[0], directions[1])
        D23 = observers[1] @ np.cross(directions[0], directions[1])
        D33 = observers[2] @ np.cross(directions[0], directions[1])

        # Parameters for calculating "r2"
        A = (-D12 * t3 / t + D22 + D32 * t1 / t) / D0
        B = (-D12 * (t**2 - t3**2) * t3 / t + D32 * (t**2 - t1**2) * t1 / t) / 6 / D0
        E = observers[1] @ directions[1]

        a = -(A**2 + 2 * A * E + np.linalg.norm(observers[1])**2)
        b = -2 * CONSTANT.GM * B * (A + E)
        c = -CONSTANT.GM**2 * B**2

        # Solve |r2|
        # Using Newton's method to Solve Equation of |r2|
        # (x^8 + a*x^6 + b*x^3+  c = 0)
        fx = np.poly1d([1, 0, a, 0, 0, b, 0, 0, c])
        dx = fx.deriv()

        r = r0
        counter = 0

        while (counter < repeated and np.abs(fx(r)) > alias):
            # print(counter)
            r -= fx(r) / dx(r)
            counter += 1

        # Find Other Constants (C1, C3)
        C1 = t3 * (1 + CONSTANT.GM / 6 / r**3 * (t**2 - t3**2)) / t
        C3 = -t1 * (1 + CONSTANT.GM / 6 / r**3 * (t**2 - t1**2)) / t

        # Find Distances from Observations to Satellite (mul_p1, mul_p2, mul_p3)
        mul_p1 = (-D11 + D21 / C1 - C3 * D31 / C1) / D0
        mul_p2 = A + CONSTANT.GM * B / r**3
        mul_p3 = (-C3 * D13 / C1 + D23 / C3 - D33) / D0

        # Find Position Vectors of Satellite (r1, r2, r3)
        r1 = observers[0] + mul_p1 * directions[0]
        r2 = observers[1] + mul_p2 * directions[1]
        r3 = observers[2] + mul_p3 * directions[2]

        # Find Velocity of Satellite at t2 (v2)
        f1 = 1 - CONSTANT.GM * t1**2 / 2 / r
        g1 = t1 - CONSTANT.GM * t1**3 / 6 / r**3
        f3 = 1 - CONSTANT.GM * t3**2 / 2 / r
        g3 = t3 - CONSTANT.GM * t3**3 / 6 / r**3

        v2 = (f1 * r3 - f3 * r1) / (f1 * g3 - f3 * g1)

        return calc_oe_from_sv(r2, v2)


"""
TISAT 1                 
1 36799U 10035E   20215.46271709  .00000102  00000-0  16380-4 0  9994
2 36799  98.1977  27.5125 0012384 172.8977 187.2418 14.91391657 545784

"""

if __name__ == "__main__":
    OE:TypeOE ={"eccentricity": 0.0012384, "inclination":98.1977, "right_ascension":27.5125, "argument_of_perigee":172.8977, "mean_anomaly":187.2418,"semimajor_axis":TypeOE.calSemimajorAxis(14.91391657)}
    observer:TypeLatLong = TypeLatLong(21.047198, 105.800237) # 18, Hoang Quoc Viet, Cau Giay, Ha Noi, Viet Nam
    time_point = datetime(year=2020, month=8, day=15, hour=10, tzinfo=None)  # Time in UTC
    simulation  = Simulation(OE, observer, time = time_point)
    orbitcalulator = OrbitCalculate(simulation)
    new_OE = orbitcalulator.calculate([3900, 4000, 4100])
    print(new_OE)