import pytest
import dataclasses
import numpy as np
import numpy.testing as npt
from datetime import datetime
from package.utility import calVecInPQW, calMeanAnomaly, calEccentricAnomaly, TypeOE, TypeLatLong
from package.orbital import Simulation, OrbitCalculate


@pytest.mark.skip
def test_Gauss():
    OE = TypeOE(eccentricity=0.0013752048537967394, semimajor_axis=6890.542965539954, inclination=1.696632359184506, right_ascension=2.7786596260229453, argument_of_perigee=0.24912363887968036, mean_anomaly=-0.18780826514804805)
    observer = TypeLatLong(21.047198, 105.800237)  # 18, Hoang Quoc Viet, Cau Giay, Ha Noi, Viet Nam
    sim = Simulation(OE, observer)
    orbit = OrbitCalculate(sim)
    orbit.observers = [np.array([-2936.2922, -5654.01, 0]), np.array([-389.5577, -6359.079, 0]), np.array([2224.1957, -5970.1419, 0])]
    orbit.directions = [np.array([-0.39868932, 0.91626844, -0.0387166]), np.array([-0.56171108, 0.82642542, 0.03875162]), np.array([-0.70363441, 0.70400468, 0.09631214])]
    time_points = np.array([3900, 4000, 4100], dtype='float64')
    OE_cal = orbit.calculate(time_points, shouldPrecalculate=False)
    OE_sample = TypeOE(eccentricity=0.0013752048537967394, semimajor_axis=6890.542965539954, inclination=1.696632359184506, right_ascension=2.7786596260229453, argument_of_perigee=0.24912363887968036, mean_anomaly=-0.18780826514804805)
    npt.assert_allclose(dataclasses.astuple(OE_cal), dataclasses.astuple(OE_sample))


@pytest.mark.skip
def test_calVecInPQW():
    OE = TypeOE(eccentricity=0.0013752048537967394, semimajor_axis=6890.542965539954, inclination=1.696632359184506, right_ascension=2.7786596260229453, argument_of_perigee=0.24912363887968036, mean_anomaly=-0.18780826514804805)
    mean_anomaly = calMeanAnomaly(OE.mean_anomaly, OE.semimajor_axis, dt=-100)
    eccentric_anomaly = calEccentricAnomaly(mean_anomaly, OE.eccentricity, 1)
    cal_PQW = calVecInPQW(OE.semimajor_axis, OE.eccentricity, eccentric_anomaly)
    sample_PQW = np.array([6577.38740598, -2023.0685358, 0])
    npt.assert_almost_equal(sample_PQW, cal_PQW, decimal=6)
