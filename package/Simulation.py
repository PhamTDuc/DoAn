import numpy as np
from typing import Tuple
from datetime import datetime
from .utility import (TypeOE, TypeLadLong, TypeXYZ, calMeanAnomaly, calVecInPQW, toECIfromLadLong, getRie)


class Simulation(object):
    def __init__(self, OE: TypeOE, observer: TypeLadLong, E0: float = 1, time: datetime = datetime.now()) -> None:
        self.OE = OE
        self.eccentric_anomaly = calEccentricAnomaly(OE['mean_anomaly'], OE['eccentricity'], E0)
        self.observer = observer  # NameTupled(ladtitude, longtitude) of Observer
        self.time = time  # Datetime is used in calculate Rie(Rotation Frame from ECEF to ECI)

    def __update(self, dt: float) -> None:
        """Recalculate MeanAbnomaly  and EccentricAnomaly after dt
        """
        OE['mean_anomaly'] = calMeanAnomaly(self.OE['mean_anomaly'], self.OE['semimajor_axis'], dt)
        self.eccentric_anomaly = calEccentricAnomaly(OE['mean_anomaly'], OE['eccentricity'], self.eccentric_anomaly)

    def __calPosPQW(dt: float) -> np.array:
        """Calculate Position Vector of Satellite in ECIF Frame
        """
        __update(dt)
        return calVecInPQW(self.OE['semimajor_axis'], self.OE['eccentricity'], self.eccentric_anomaly)

    def _getObserverECI(dt: float) -> TypeXYZ:
        return toECIfromLadLong(self.observer.ladtitude, self.observer.longtitude, time=self.time + datetime.second(dt))

    def _getPosECI(dt: float) -> TypeXYZ:
        pos_pqw = __calPosPQW(dt)
        toECI = getMatECItoPQW(i=self.OE['inclination'], omega=self.OE['argument_of_perigee'], sigma=self.OE['right_ascension'])
        coord = toECI.T @ pos_pqw
        return TypeXYZ(*coord)

    def getRie(dt: float) -> np.array:
        """Calculate Rie from self.time + dt(second)
        """
        return getRie(self.time + datetime.second(dt))

    def getAllCoords(dt: float) -> Tuple[TypeXYZ, TypeXYZ]:
        pos = _getPosECI(dt)
        observer = _getObserverECI(dt)
        return observer, pos
