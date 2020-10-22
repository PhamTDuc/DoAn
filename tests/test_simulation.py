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

@pytest.mark.xfail
@pytest.mark.parametrize("dt, true_result",[(3900, np.array([-2936.2922, -5654.01, 0])),\
(4000,np.array([-389.5577, -6359.079, 0])),\
(4100,np.array([2224.1957, -5970.1419, 0]))])
def test_calulateObservers(getSimulation, dt:float, true_result):

    sim = getSimulation
    observer = sim.getAllCoords(dt)[0]
    npt.assert_almost_equal(observer, true_result, decimal=2)
