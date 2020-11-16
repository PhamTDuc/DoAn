import pytest
import numpy as np
import numpy.testing as npt
from datetime import datetime, timedelta
from package.utility import TypeOE, getGMST


def test_SemimajorAxisfromMeanMotion():
    OE = TypeOE(mean_motion=15.5918272)
    npt.assert_allclose(OE.semimajor_axis, 6768.0201)


@pytest.mark.skip
@pytest.mark.parametrize("dt, expected",[(3900, 136.18028935), (4000, 136.598096813),(4100, 137.0159042756)])
def test_GMST(dt, expected):
	time = datetime(year = 2020, month = 4, day= 15, hour=18, minute=22, second=4)
	npt.assert_almost_equal(getGMST(time+timedelta(seconds=dt)), expected, decimal=2)