from __future__ import annotations

from importlib.resources import files
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator


class SecondaryLossSurrogate:
    """304-strand secondary-subcable full-pulse loss surrogate.

    The database represents a fixed 3/1/3 ms trapezoidal pulse. Directional
    dependence is stored as the quadratic form

        Q(theta)=Qxx cos^2(theta)+Qyy sin^2(theta)+2 Qxy sin(theta)cos(theta).

    Below the first positive field node, Q/B^2 is held constant to preserve the
    correct low-field scaling Q~B^2.
    """

    def __init__(self, csv_path: str | Path | None = None):
        if csv_path is None:
            csv_path = files("apexloss.data").joinpath("secondary304_surrogate_grid.csv")
        self.df = pd.read_csv(csv_path)
        self.T_grid = np.sort(self.df["T_K"].unique())
        self.B_grid = np.sort(self.df["Bpk_T"].unique())

    def _component_at_T(self, col, Bq, Tnode):
        d = self.df[np.isclose(self.df["T_K"], Tnode)].sort_values("Bpk_T")
        b = d["Bpk_T"].to_numpy(float)
        q = d[col].to_numpy(float)
        m = b > 0
        bp = b[m]
        s = q[m] / bp**2
        f = PchipInterpolator(bp, s, extrapolate=True)
        Bq = np.asarray(Bq, float)
        out = np.zeros_like(Bq)
        pos = Bq > 0
        bc = np.clip(Bq[pos], bp[0], bp[-1])
        sv = f(bc)
        sv = np.where(Bq[pos] < bp[0], s[0], sv)
        sv = np.where(Bq[pos] > bp[-1], s[-1], sv)
        out[pos] = Bq[pos] ** 2 * sv
        return out

    def _component(self, col, Bq, Tq):
        Bq, Tq = np.broadcast_arrays(np.asarray(Bq, float), np.asarray(Tq, float))
        vals = np.stack([self._component_at_T(col, Bq, T) for T in self.T_grid], axis=0)
        out = np.empty(Bq.size)
        Tf = Tq.ravel()
        flat = vals.reshape(len(self.T_grid), -1)
        for i in range(Bq.size):
            out[i] = np.interp(Tf[i], self.T_grid, flat[:, i])
        return out.reshape(Bq.shape)

    def q_secondary_J_per_m(self, Br_peak_T, Bz_peak_T, T_K=4.2):
        Br, Bz, T = np.broadcast_arrays(
            np.asarray(Br_peak_T, float), np.asarray(Bz_peak_T, float), np.asarray(T_K, float)
        )
        B = np.hypot(Br, Bz)
        th = np.arctan2(Bz, Br)
        c = np.cos(th)
        s = np.sin(th)
        qxx = self._component("Qxx_Jpm", B, T)
        qyy = self._component("Qyy_Jpm", B, T)
        qxy = self._component("Qxy_Jpm", B, T)
        q = qxx * c**2 + qyy * s**2 + 2 * qxy * s * c
        return float(q) if q.ndim == 0 else q

    def q_tertiary_J_per_m(self, Br_peak_T, Bz_peak_T, T_K=4.2, n_secondary=8):
        return n_secondary * self.q_secondary_J_per_m(Br_peak_T, Bz_peak_T, T_K)
