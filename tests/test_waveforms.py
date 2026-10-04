import numpy as np
from apexloss.waveforms import TrapezoidPulse


def test_trapezoid_breakpoints():
    p = TrapezoidPulse()
    assert p.amplitude(0.0) == 0.0
    assert np.isclose(p.amplitude(1.5e-3), 0.5)
    assert np.isclose(p.amplitude(3e-3), 1.0)
    assert np.isclose(p.amplitude(4e-3), 1.0)
    assert np.isclose(p.amplitude(5.5e-3), 0.5)
    assert np.isclose(p.amplitude(7e-3), 0.0)
