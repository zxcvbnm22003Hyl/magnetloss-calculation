"""Run the compact release-validation suite without pytest."""
from __future__ import annotations

import json
import math

from apexloss import MagnetGeometry, FiniteTurnBiotSavart, SingleStrandDiffusion, SecondaryLossSurrogate
from apexloss.materials import LegacyKohlerRRR
from apexloss.whole_magnet import fixed_temperature_loss


def main():
    strand = SingleStrandDiffusion(0.2e-3, 160)
    r = strand.simulate_trapezoid(
        B_peak_T=20.619192, dt_s=1e-6, tail_s=23e-3, resistivity=LegacyKohlerRRR()
    )
    assert math.isclose(r.Q_pulse_J_per_m, 0.5810907987, rel_tol=2e-5)

    field = FiniteTurnBiotSavart(MagnetGeometry(), 47.29e3, 20)
    assert math.isclose(field.center_field_T(), 20.1952285370, rel_tol=2e-6)

    _, q = fixed_temperature_loss(field, SecondaryLossSurrogate(), 4.2, 2)
    assert math.isclose(q, 73335.74799, rel_tol=5e-4)

    print(json.dumps({
        "strand_Qpulse_Jpm": r.Q_pulse_J_per_m,
        "strand_Qtotal_Jpm": r.Q_total_J_per_m,
        "B0_T": field.center_field_T(),
        "whole_magnet_fixed4p2K_Qeddy_J": q,
        "status": "PASS"
    }, indent=2))


if __name__ == "__main__":
    main()
