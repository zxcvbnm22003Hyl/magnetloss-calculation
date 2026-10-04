import numpy as np
from apexloss import SecondaryLossSurrogate


def test_surrogate_reference_point():
    s = SecondaryLossSurrogate()
    q = s.q_secondary_J_per_m(20.5, 0.0, 4.2)
    assert np.isclose(q, 178.666134, rtol=3e-3)


def test_low_field_quadratic_scaling():
    s = SecondaryLossSurrogate()
    q1 = s.q_secondary_J_per_m(0.5, 0.0, 4.2)
    q2 = s.q_secondary_J_per_m(1.0, 0.0, 4.2)
    assert np.isclose(q2 / q1, 4.0, rtol=1e-12)
