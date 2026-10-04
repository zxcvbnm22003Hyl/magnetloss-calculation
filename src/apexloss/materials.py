from __future__ import annotations

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import cumulative_trapezoid


class ResistivityTable:
    """Bilinear rho(T,B) interpolator for the bundled high-purity Al model.

    The repository stores only the compact zero-field rho0(T) curve. At
    initialization the original 0--20.5 T, 0.5 T-spaced grid is reconstructed
    with the same Kohler-type magnetoresistance law used to build the historical
    full table.
    """

    rho293_ohm_m: float = 2.94978e-8

    def __init__(self, csv_path: str | Path | None = None, clip: bool = True):
        if csv_path is None:
            csv_path = files("apexloss.data").joinpath("rho0_al_r3000_T.csv")
        df = pd.read_csv(csv_path)
        if {"B_T", "rho_TB_Ohm_m"}.issubset(df.columns):
            self.T_grid = np.sort(df["T_K"].unique())
            self.B_grid = np.sort(df["B_T"].unique())
            table = (
                df.pivot(index="T_K", columns="B_T", values="rho_TB_Ohm_m")
                .loc[self.T_grid, self.B_grid]
                .to_numpy(float)
            )
        else:
            if not {"T_K", "rho0_Ohm_m"}.issubset(df.columns):
                raise ValueError(
                    "Material CSV must contain either (T_K, B_T, rho_TB_Ohm_m) "
                    "or compact (T_K, rho0_Ohm_m) columns."
                )
            compact = df[["T_K", "rho0_Ohm_m"]].drop_duplicates().sort_values("T_K")
            self.T_grid = compact["T_K"].to_numpy(float)
            rho0 = compact["rho0_Ohm_m"].to_numpy(float)
            self.B_grid = np.arange(0.0, 20.5 + 0.25, 0.5)
            K_rho = self.rho293_ohm_m / rho0[:, None]
            H = 0.01 * self.B_grid[None, :] * K_rho
            MR = H**2 * (1.0 + 0.00177 * H) / (1.8 + 1.6 * H + 0.53 * H**2)
            table = rho0[:, None] * (1.0 + MR)

        self._interp = RegularGridInterpolator(
            (self.T_grid, self.B_grid), table, bounds_error=False, fill_value=None
        )
        self.clip = clip

    def __call__(self, T_K, B_T):
        T = np.asarray(T_K, dtype=float)
        B = np.asarray(B_T, dtype=float)
        shape = np.broadcast_shapes(T.shape, B.shape)
        T = np.broadcast_to(T, shape).ravel()
        B = np.broadcast_to(B, shape).ravel()
        if self.clip:
            T = np.clip(T, self.T_grid[0], self.T_grid[-1])
            B = np.clip(B, self.B_grid[0], self.B_grid[-1])
        y = self._interp(np.column_stack([T, B])).reshape(shape)
        return float(y) if y.ndim == 0 else y


@dataclass(frozen=True)
class LegacyKohlerRRR:
    """Legacy 4.2 K magnetoresistance fit used by the original strand model."""

    rho0_ohm_m: float = 9.832581476e-12
    rho293_ohm_m: float = 2.94978e-8

    def __call__(self, B_T):
        B = np.asarray(B_T, dtype=float)
        K_rho = self.rho293_ohm_m / self.rho0_ohm_m
        H = 0.01 * np.abs(B) * K_rho
        MR = H**2 * (1 + 0.00177 * H) / (1.8 + 1.6 * H + 0.53 * H**2)
        out = self.rho0_ohm_m * (1 + MR)
        return float(out) if out.ndim == 0 else out


def cp_aluminum_proxy(T_K):
    """Al-3003-F proxy heat capacity polynomial retained for model continuity."""
    T = np.maximum(np.asarray(T_K, dtype=float), 4.0)
    x = np.log10(T)
    y = (
        46.6467
        - 314.292 * x
        + 866.662 * x**2
        - 1298.3 * x**3
        + 1162.27 * x**4
        - 637.795 * x**5
        + 210.351 * x**6
        - 38.3094 * x**7
        + 2.96344 * x**8
    )
    return 10**y


class EnthalpyModel:
    """Numerical h(T) and inverse T(h) map for adiabatic pulse calculations."""

    def __init__(self, Tmin_K=4.2, Tmax_K=293.15, n=20000):
        self.T_grid = np.linspace(Tmin_K, Tmax_K, n)
        self.h_grid = np.concatenate(
            [[0.0], cumulative_trapezoid(cp_aluminum_proxy(self.T_grid), self.T_grid)]
        )

    def temperature_from_specific_enthalpy(self, h_J_kg):
        return np.interp(h_J_kg, self.h_grid, self.T_grid)
