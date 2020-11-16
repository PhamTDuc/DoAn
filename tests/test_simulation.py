import pytest
from datetime import datetime
from package.utility import TypeOE, TypeLatLong
from package.orbital import Simulation
import numpy as np
import numpy.testing as npt

@pytest.fixture
def getSimulation():
    OE = TypeOE(semimajor_axis=6876.6644, eccentricity=0.0020122, inclination=1.6971285, \
    right_ascension=2.7791209, argument_of_perigee=1.6541096, mean_anomaly=0.2653564)
    observer = TypeLatLong(0, 0)
    time = datetime(year=2020, month=4, day=15, hour=18, minute=22, second=4)

    return Simulation(OE, observer, time=time)

@pytest.mark.skip
@pytest.mark.parametrize("dt, true_result",[(3900, np.array([-4596.8172, 4411.2258, 0])),\
(4000,np.array([-4628.8618, 4377.5883, 0])),\
(4100,np.array([-4660.6604, 4343.7180, 0]))])
def test_calulateObservers(getSimulation, dt:float, true_result):

    sim = getSimulation
    observer = sim.getAllCoords(dt)[0]
    npt.assert_almost_equal(observer, true_result, decimal=2)


@pytest.mark.parametrize("dt, true_result",[(3900, np.array([-1842.0394, -2015.6368, -340.1330])),\
(4000,np.array([-1770.6852, -1894.3849, -414.6002])),\
(4100,np.array([-1621.2242, -1813.3035, 1164.2916]))])
def test_calculateDirectionVectors(getSimulation, dt:float, true_result):
    sim = getSimulation
    direction = sim.getAllCoords(dt)[1]
    npt.assert_almost_equal(direction, true_result, decimal=2)
