from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.special import ellipk, ellipe
from numpy.polynomial.legendre import leggauss

from .constants import MU0
from .geometry import MagnetGeometry


def circular_loop_field(radius_m, source_z_m, current_A, r_m, z_m):
    """Magnetic field of one or many circular filamentary loops.

    Parameters are broadcast-compatible except the observation point, which is scalar.
    Returns Br, Bz arrays corresponding to source loops.
    """
    a = np.asarray(radius_m, dtype=float)
    z0 = np.asarray(source_z_m, dtype=float)
    I = np.asarray(current_A, dtype=float)
    r = float(r_m)
    dz = float(z_m) - z0

    if abs(r) < 1e-14:
        Br = np.zeros_like(a)
        Bz = MU0 * I * a**2 / (2 * (a**2 + dz**2) ** 1.5)
        return Br, Bz

    den = (a + r) ** 2 + dz**2
    root = np.sqrt(den)
    delta = (a - r) ** 2 + dz**2
    m = np.clip(4 * a * r / den, 0.0, 1.0 - 1e-14)
    K = ellipk(m)
    E = ellipe(m)
    common = MU0 * I / (2 * np.pi * root)
    Br = common * (dz / r) * (-K + (a**2 + r**2 + dz**2) / delta * E)
    Bz = common * (K + (a**2 - r**2 - dz**2) / delta * E)
    return Br, Bz


@dataclass
class FiniteTurnBiotSavart:
    """Axisymmetric finite rectangular-turn Biot-Savart model.

    Each macro turn is integrated as a continuum of circular loops by tensor-product
    Gauss-Legendre quadrature. Use an even source order to avoid a source point exactly
    at a turn center when evaluating self-field there.
    """

    geometry: MagnetGeometry = MagnetGeometry()
    current_A: float = 47.29e3
    source_order: int = 20

    def __post_init__(self):
        if self.source_order % 2:
            raise ValueError("source_order should be even to avoid center self-singular sampling")
        self._build_sources()

    def _build_sources(self):
        x, w = leggauss(self.source_order)
        g = self.geometry
        Jphi = self.current_A / (g.cable_radial_m * g.cable_axial_m)
        aa, zz, ii = [], [], []
        for zc in g.axial_centers():
            for rc in g.radial_centers():
                aq = rc + 0.5 * g.cable_radial_m * x
                zq = zc + 0.5 * g.cable_axial_m * x
                wa = 0.5 * g.cable_radial_m * w
                wz = 0.5 * g.cable_axial_m * w
                A, Z = np.meshgrid(aq, zq, indexing="ij")
                WA, WZ = np.meshgrid(wa, wz, indexing="ij")
                aa.append(A.ravel())
                zz.append(Z.ravel())
                ii.append((Jphi * WA * WZ).ravel())
        self._a = np.concatenate(aa)
        self._z = np.concatenate(zz)
        self._I = np.concatenate(ii)

    def field(self, r_m: float, z_m: float) -> tuple[float, float]:
        Br, Bz = circular_loop_field(self._a, self._z, self._I, r_m, z_m)
        return float(np.sum(Br)), float(np.sum(Bz))

    def center_field_T(self) -> float:
        return self.field(0.0, 0.0)[1]

    def turn_center_fields(self) -> pd.DataFrame:
        table = self.geometry.turn_table()
        br, bz = [], []
        for row in table.itertuples(index=False):
            x, y = self.field(row.r_m, row.z_m)
            br.append(x)
            bz.append(y)
        table["Br_pk_T"] = br
        table["Bz_pk_T"] = bz
        table["Bmag_pk_T"] = np.hypot(table["Br_pk_T"], table["Bz_pk_T"])
        return table

    def turn_gauss_points(self, order: int = 2) -> pd.DataFrame:
        """Return field at area-normalized Gauss points in every turn cross section."""
        x, w = leggauss(order)
        g = self.geometry
        rows = []
        for iz, zc in enumerate(g.axial_centers()):
            for ir, rc in enumerate(g.radial_centers()):
                turn_id = iz * g.n_radial + ir + 1
                for i, xi in enumerate(x):
                    for j, eta in enumerate(x):
                        r = rc + 0.5 * g.cable_radial_m * xi
                        z = zc + 0.5 * g.cable_axial_m * eta
                        wt = w[i] * w[j] / 4.0
                        Br, Bz = self.field(r, z)
                        rows.append(
                            {
                                "turn_id": turn_id,
                                "sp_id": iz + 1,
                                "dp_id": iz // 2 + 1,
                                "radial_layer": ir + 1,
                                "r_center_m": rc,
                                "z_center_m": zc,
                                "r_sample_m": r,
                                "z_sample_m": z,
                                "weight": wt,
                                "Br_pk_T": Br,
                                "Bz_pk_T": Bz,
                                "Bmag_pk_T": float(np.hypot(Br, Bz)),
                                "turn_length_m": 2 * np.pi * rc,
                            }
                        )
        return pd.DataFrame(rows)
