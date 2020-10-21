import numpy as np
from datetime import datetime, timedelta
from typing import Tuple
from .utility import (CONSTANT, TypeOE, TypeLatLong, TypeXYZ, normalize, calMeanAnomaly, calEccentricAnomaly, calVecInPQW, toECIfromLatLong, getRie, getMatECItoPQW)
from .orbital_elements import calc_oe_from_sv
import dataclasses


class Simulation(object):
    def __init__(self, OE: TypeOE, observer: TypeLatLong, E0: float = 1, time: datetime = datetime.now()) -> None:
        self._OE = dataclasses.replace(OE)
        self.eccentric_anomaly = None
        self._observer = observer
        self._time = time

    def _update(self, dt: float) -> None:
        self._OE.mean_anomaly = calMeanAnomaly(self._OE.mean_anomaly, self._OE.semimajor_axis, dt)
        self.eccentric_anomaly = calEccentricAnomaly(self._OE.mean_anomaly, self._OE.eccentricity, alias=1)

    def _getObserverECI(self, dt: float) -> np.array:
        return toECIfromLatLong(self._observer, time=self._time + timedelta(seconds=dt))

    def _getPosInECI(self) -> np.array:
        pos_pqw = calVecInPQW(self._OE.semimajor_axis, self._OE.eccentricity, self.eccentric_anomaly)
        toECI = getMatECItoPQW(i=self._OE.inclination, omega=self._OE.argument_of_perigee, sigma=self._OE.right_ascension)
        return toECI.T @ pos_pqw

    def getAllCoords(self, dt: float) -> Tuple[np.array, np.array]:
        self._update(dt)
        pos = self._getPosInECI()
        observer = self._getObserverECI(dt)
        return observer, pos


class OrbitCalculate(object):

    def __init__(self, sim: Simulation) -> None:
        self.sim = sim
        self.observers = []
        self.directions = []

    def preCalculate(self, time_points: Tuple[float, float, float]):
        self.observers = []
        self.directions = []

        for dt in time_points:
            observer, pos = self.sim.getAllCoords(dt)
            direction = normalize(pos - observer)
            self.observers.append(observer)
            self.directions.append(direction)

    def calculate(self, time_points: Tuple[float, float, float], r0: float=600, repeated: int=10000, alias: float = 0.1, shouldPrecalculate: bool=True) -> TypeOE:

        if (repeated < 0 or alias < 0):
            raise ValueError("Repeated must be int. Repeated and Alias must be positive.")

        if(shouldPrecalculate):
            preCalculate(time_points)

        assert len(self.observers) == 3, "Obeservers must be array of 3 Vec3"
        assert len(self.directions) == 3, "Obeservers must be array of 3 Vec3"
        t1 = time_points[0] - time_points[1]
        t3 = time_points[2] - time_points[1]
        t = t3 - t1

        D0 = self.directions[0] @ np.cross(self.directions[1], self.directions[2])
        D11 = self.observers[0] @ np.cross(self.directions[1], self.directions[2])
        D21 = self.observers[1] @ np.cross(self.directions[1], self.directions[2])
        D31 = self.observers[2] @ np.cross(self.directions[1], self.directions[2])

        D12 = self.observers[0] @ np.cross(self.directions[0], self.directions[2])
        D22 = self.observers[1] @ np.cross(self.directions[0], self.directions[2])
        D32 = self.observers[2] @ np.cross(self.directions[0], self.directions[2])

        D13 = self.observers[0] @ np.cross(self.directions[0], self.directions[1])
        D23 = self.observers[1] @ np.cross(self.directions[0], self.directions[1])
        D33 = self.observers[2] @ np.cross(self.directions[0], self.directions[1])

        # Parameters for calculating "r2"
        A = (-D12 * t3 / t + D22 + D32 * t1 / t) / D0
        B = (-D12 * (t**2 - t3**2) * t3 / t + D32 * (t**2 - t1**2) * t1 / t) / 6 / D0
        E = self.observers[1] @ self.directions[1]

        a = -(A**2 + 2 * A * E + np.linalg.norm(self.observers[1])**2)
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
        r1 = self.observers[0] + mul_p1 * self.directions[0]
        r2 = self.observers[1] + mul_p2 * self.directions[1]
        r3 = self.observers[2] + mul_p3 * self.directions[2]

        # Find Velocity of Satellite at t2 (v2)
        f1 = 1 - CONSTANT.GM * t1**2 / 2 / r
        g1 = t1 - CONSTANT.GM * t1**3 / 6 / r**3
        f3 = 1 - CONSTANT.GM * t3**2 / 2 / r
        g3 = t3 - CONSTANT.GM * t3**3 / 6 / r**3

        v2 = (f1 * r3 - f3 * r1) / (f1 * g3 - f3 * g1)

        return calc_oe_from_sv(r2, v2)


if __name__ == "__main__":
    OE = TypeOE(eccentricity=0.0012384, inclination=98.1977, right_ascension=27.5125, argument_of_perigee=172.8977, mean_anomaly=187.2418, mean_motion=14.91391657)
    print(OE)
