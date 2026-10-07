import numpy as np
import pytest

from apexloss.complete_penetration import CompletePenetrationStrand
from apexloss.materials import LegacyKohlerRRR
from apexloss.waveforms import TrapezoidPulse


def test_complete_penetration_constant_rho_matches_closed_form():
    model = CompletePenetrationStrand(diameter_m=0.2e-3)
    pulse = TrapezoidPulse()
    rho = 2.5e-9
    Bpk = 12.0

    q_exact = model.closed_form_constant_rho(Bpk, rho, pulse)
    coeff = np.pi * model.radius_m**4 / (4.0 * rho)
    q_manual = coeff * Bpk**2 * (
        1.0 / pulse.rise_s + 1.0 / pulse.decay_s
    )

    assert q_exact == pytest.approx(q_manual, rel=1e-14)


def test_complete_penetration_legacy_regression():
    model = CompletePenetrationStrand(diameter_m=0.2e-3)
    r = model.simulate_trapezoid(
        B_peak_T=20.619192,
        dt_s=1e-6,
        resistivity=LegacyKohlerRRR(),
    )

    assert r.Q_pulse_J_per_m == pytest.approx(0.6083367961, rel=2e-8)
    assert r.Q_tail_J_per_m == 0.0
    assert r.Q_total_J_per_m == pytest.approx(r.Q_pulse_J_per_m)


def test_complete_penetration_flat_top_power_is_zero():
    model = CompletePenetrationStrand(diameter_m=0.2e-3)
    p = model.power(rho_ohm_m=1e-9, dBdt_T_s=0.0)
    assert float(p) == 0.0
