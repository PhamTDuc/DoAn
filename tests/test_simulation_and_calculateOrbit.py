import pytest
from datetime import datetime
from package.utility import TypeOE, TypeLatLong, normalize
from package.orbital import Simulation, OrbitCalculate
import numpy as np
import numpy.testing as npt


@pytest.mark.skip
@pytest.fixture
def getSimulation():
    OE = TypeOE(semimajor_axis=6876.6644, eccentricity=0.0020122, inclination=1.6971285,
                right_ascension=2.7791209, argument_of_perigee=1.6541096, mean_anomaly=0.2653564)
    observer = TypeLatLong(0, 0)
    time = datetime(year=2020, month=4, day=15, hour=18, minute=22, second=4)

    return Simulation(OE, observer, time=time)


@pytest.mark.skip
@pytest.mark.parametrize("dt, true_result", [(3900, np.array([-4596.8172, 4411.2258, 0])),
                                             (4000, np.array([-4628.8618, 4377.5883, 0])),
                                             (4100, np.array([-4660.6604, 4343.7180, 0]))])
def test_calulateObservers(getSimulation, dt: float, true_result):

    sim = getSimulation
    observer = sim.getAllCoords(dt)[0]
    npt.assert_almost_equal(observer, true_result, decimal=2)


@pytest.mark.skip
@pytest.mark.parametrize("dt, true_result", [(3900, np.array([-1842.0394, -2015.6368, -340.1330])),
                                             (4000, np.array([-1770.6852, -1894.3849, 414.6002])),
                                             (4100, np.array([-1621.2242, -1813.3035, 1164.2916]))])
def test_calculateDirectionVectors(getSimulation, dt: float, true_result):
    sim = getSimulation
    direction = sim.getAllCoords(dt)[1]
    npt.assert_almost_equal(direction, true_result, decimal=2)


@pytest.mark.skip
@pytest.mark.parametrize("dt, true_result", [(3900, np.array([-0.6694, -0.7325, -0.1236])),
                                             (4000, np.array([-0.6743, -0.7214, 0.1579])),
                                             (4100, np.array([-0.6027, -0.6704, 0.4328]))])
def test_DirectionNormalize(getSimulation, dt: float, true_result):
    sim = getSimulation
    direction = normalize(sim.getAllCoords(dt)[1])
    npt.assert_almost_equal(direction, true_result, decimal=2)


# @pytest.mark.skip
def test_Gauss(getSimulation):
    sim = getSimulation
    orbit = OrbitCalculate(sim)
    # orbit.observers = [np.array([-4596.8172, 4411.2258, 0]), np.array([-4628.8618, 4377.5883, 0]), np.array([-4660.6604, 4343.7180, 0])]
    # orbit.directions = [np.array([-0.6694, -0.7325, -0.1236]), np.array([-0.6743, -0.7214, 0.1579]), np.array([-0.6027, -0.6704, 0.4328])]
    Calculated_OE = orbit.calculate([3900, 4000, 4100], shouldPrecalculate=True)
    Expedted_OE = TypeOE(eccentricity=0.009431418961541798, semimajor_axis=6809.778031000419, inclination=1.696733706840077, right_ascension=2.7775731301469375, argument_of_perigee=3.006356982240637, mean_anomaly=-2.942070455156746)
    assert(Calculated_OE == Expedted_OE)
