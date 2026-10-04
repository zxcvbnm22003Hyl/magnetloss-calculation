import pytest
from apexloss.strand_diffusion import SingleStrandDiffusion
from apexloss.materials import LegacyKohlerRRR


@pytest.mark.slow
def test_legacy_single_strand_regression():
    model = SingleStrandDiffusion(diameter_m=0.2e-3, radial_cells=160)
    r = model.simulate_trapezoid(
        B_peak_T=20.619192,
        dt_s=1e-6,
        tail_s=23e-3,
        resistivity=LegacyKohlerRRR(),
    )
    assert r.Q_pulse_J_per_m == pytest.approx(0.5810907987, rel=2e-5)
    assert r.Q_tail_J_per_m == pytest.approx(0.0072395840, rel=2e-5)
    assert r.Q_total_J_per_m == pytest.approx(0.5883303827, rel=2e-5)
    assert r.Q_complete_penetration_J_per_m == pytest.approx(0.6083368195, rel=2e-5)
