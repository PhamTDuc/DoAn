import unittest
from datetime import datetime
import numpy as np
import numpy.testing as npt
from package.utility import calVecInPQW, calMeanAnomaly, calEccentricAnomaly, TypeOE, TypeLatLong
from package.orbital import Simulation


class TestUtility(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        print(f"##----Starting--{cls.__name__}")
        cls.OE = TypeOE(eccentricity=0.0013752048537967394, semimajor_axis=6890.542965539954, inclination=1.696632359184506, right_ascension=2.7786596260229453, argument_of_perigee=0.24912363887968036, mean_anomaly=-0.18780826514804805)
        cls.sample_PQW = np.array([6577.38740598, -2023.0685358, 0])

    @classmethod
    def tearDownClass(cls):
        print(f"##----End---{cls.__name__}")

    def test_calVecInPQW(self):
        mean_anomaly = calMeanAnomaly(type(self).OE.mean_anomaly, type(self).OE.semimajor_axis, dt=-100)
        eccentric_anomaly = calEccentricAnomaly(mean_anomaly, type(self).OE.eccentricity, 1)
        vec_PQW = calVecInPQW(type(self).OE.semimajor_axis, type(self).OE.eccentricity, eccentric_anomaly)
        npt.assert_almost_equal(type(self).sample_PQW, vec_PQW, decimal=6)

    def test_getPosECI(self):
        observer: TypeLatLong = TypeLatLong(21.047198, 105.800237)  # 18, Hoang Quoc Viet, Cau Giay, Ha Noi, Viet Nam
        time_point = datetime(year=2020, month=8, day=15, hour=10, tzinfo=None)  # Time in UTC
        simulation = Simulation(type(self).OE, observer, time=time_point)
        simulation._update(3900)
        ret = simulation._getPosInECI()
        print(ret)


class TestFunction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        print(f"##----Starting--{cls.__name}")
        OE = TypeOE(eccentricity=0.0012384, inclination=98.1977, right_ascension=27.5125, argument_of_perigee=172.8977, mean_anomaly=187.2418, mean_motion=14.91391657)
        observer = TypeLatLong(21.047198, 105.800237)  # 18, Hoang Quoc Viet, Cau Giay, Ha Noi, Viet Nam
        sim = Simulation(OE, observer)
        cls.orbit = OrbitCalculate(sim)

    @classmethod
    def tearDownClass(cls):
        print(f"##----End---{cls.__name__}")

    def test_Gauss(self):
        time_points = np.array([3900, 4000, 4100], dtype='float64')
        new_OE = type(self).orbit.calculate(time_points)
        npt.assert_approx_allclose(np.array([new_OE.eccentricity, new_OE.semimajor_axis, new_OE.inclination, new_OE.right_ascension, new_OE.argument_of_perigee, new_OE.mean_anomaly]),
            np.array([0.0013752048537967394, 6890.542965539954, 1.696632359184506, 2.7786596260229453, 0.24912363887968036, -0.18780826514804805])

if __name__ == "__main__":
    unittest.main()
