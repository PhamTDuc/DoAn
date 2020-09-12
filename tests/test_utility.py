import unittest
import numpy as np
import numpy.testing as npt
from package.utility import calVecInPQW, calMeanAnomaly, calEccentricAnomaly, TypeOE


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


if __name__ == "__main__":
    unittest.main()
