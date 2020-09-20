import pytest
import numpy as np
import numpy.testing as npt
from datetime import datetime
from package.utility import TypeOE


def test_SemimajorAxisfromMeanMotion():
    OE = TypeOE(mean_motion=15.5918272)
    npt.assert_allclose(OE.semimajor_axis, 6768.0201)
