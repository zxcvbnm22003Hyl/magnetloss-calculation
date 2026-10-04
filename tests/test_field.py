import numpy as np
from apexloss import MagnetGeometry, FiniteTurnBiotSavart


def test_phase1_center_field_regression():
    f = FiniteTurnBiotSavart(MagnetGeometry(), current_A=47.29e3, source_order=12)
    assert np.isclose(f.center_field_T(), 20.1952285, rtol=2e-6)
