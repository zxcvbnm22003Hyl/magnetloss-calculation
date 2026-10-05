import pytest

from apexloss.materials import LegacyKohlerRRR
from apexloss.modal_diffusion import RobinBesselStrand


def test_robin_bessel_legacy_regression():
    model = RobinBesselStrand(diameter_m=0.2e-3, n_modes=20)
    r = model.simulate_trapezoid(
        B_peak_T=20.619192,
        dt_s=10e-6,
        tail_s=23e-3,
        resistivity=LegacyKohlerRRR(),
    )
    assert r.Q_pulse_J_per_m == pytest.approx(0.58107356, rel=2e-5)
    assert r.Q_tail_J_per_m == pytest.approx(0.00721641, rel=2e-5)
    assert r.Q_total_J_per_m == pytest.approx(0.58828998, rel=2e-5)


def test_complete_penetration_power_is_positive():
    p = RobinBesselStrand.complete_penetration_power_W_per_m(
        radius_m=0.1e-3,
        rho_ohm_m=1e-9,
        dBdt_T_s=1000.0,
    )
    assert p > 0.0
