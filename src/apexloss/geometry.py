from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class MagnetGeometry:
    """Axisymmetric APEX-style double-pancake solenoid geometry."""

    n_dp: int = 9
    n_radial: int = 9
    inner_radius_m: float = 125e-3
    cable_radial_m: float = 8.3e-3
    cable_axial_m: float = 16.9e-3
    radial_gap_m: float = 2.0e-3
    dp_internal_gap_m: float = 0.5e-3
    inter_dp_gap_m: float = 3.0e-3

    @property
    def n_sp(self) -> int:
        return 2 * self.n_dp

    @property
    def n_turns(self) -> int:
        return self.n_sp * self.n_radial

    @property
    def outer_radius_m(self) -> float:
        return (
            self.inner_radius_m
            + self.n_radial * self.cable_radial_m
            + (self.n_radial - 1) * self.radial_gap_m
        )

    @property
    def active_length_m(self) -> float:
        return (
            self.n_dp * (2 * self.cable_axial_m + self.dp_internal_gap_m)
            + (self.n_dp - 1) * self.inter_dp_gap_m
        )

    def radial_centers(self) -> np.ndarray:
        pitch = self.cable_radial_m + self.radial_gap_m
        return self.inner_radius_m + self.cable_radial_m / 2 + np.arange(self.n_radial) * pitch

    def axial_centers(self) -> np.ndarray:
        z = []
        zcur = -self.active_length_m / 2
        for i in range(self.n_dp):
            z.append(zcur + self.cable_axial_m / 2)
            zcur += self.cable_axial_m + self.dp_internal_gap_m
            z.append(zcur + self.cable_axial_m / 2)
            zcur += self.cable_axial_m
            if i < self.n_dp - 1:
                zcur += self.inter_dp_gap_m
        return np.asarray(z)

    def turn_table(self) -> pd.DataFrame:
        rows = []
        r = self.radial_centers()
        z = self.axial_centers()
        for iz, zc in enumerate(z):
            for ir, rc in enumerate(r):
                rows.append(
                    {
                        "turn_id": iz * self.n_radial + ir + 1,
                        "sp_id": iz + 1,
                        "dp_id": iz // 2 + 1,
                        "radial_layer": ir + 1,
                        "r_m": rc,
                        "z_m": zc,
                        "turn_length_m": 2 * np.pi * rc,
                    }
                )
        return pd.DataFrame(rows)
