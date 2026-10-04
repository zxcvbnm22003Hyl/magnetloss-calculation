from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from .field import FiniteTurnBiotSavart
from .materials import ResistivityTable, EnthalpyModel
from .surrogate import SecondaryLossSurrogate
from .waveforms import TrapezoidPulse


@dataclass
class MagnetLossSummary:
    Q_eddy_J: float
    Q_transport_J: float
    Q_total_J: float
    T_mean_end_K: float | None = None
    T_max_end_K: float | None = None


def fixed_temperature_loss(
    field_model: FiniteTurnBiotSavart,
    surrogate: SecondaryLossSurrogate,
    temperature_K=4.2,
    target_order=2,
    n_secondary=8,
):
    s = field_model.turn_gauss_points(target_order)
    q2 = surrogate.q_secondary_J_per_m(
        s["Br_pk_T"].to_numpy(), s["Bz_pk_T"].to_numpy(), temperature_K
    )
    s["q_secondary_Jpm"] = q2
    s["q_tertiary_Jpm"] = n_secondary * q2
    turn_rows = []
    for tid, g in s.groupby("turn_id"):
        q = float(np.sum(g["weight"] * g["q_tertiary_Jpm"]))
        L = float(g["turn_length_m"].iloc[0])
        row = g.iloc[0]
        turn_rows.append(
            {
                "turn_id": int(tid),
                "sp_id": int(row.sp_id),
                "dp_id": int(row.dp_id),
                "radial_layer": int(row.radial_layer),
                "r_m": float(row.r_center_m),
                "z_m": float(row.z_center_m),
                "q_eddy_Jpm": q,
                "Q_eddy_J": q * L,
            }
        )
    turns = pd.DataFrame(turn_rows)
    return turns, float(turns["Q_eddy_J"].sum())


class QuasiStaticThermalMapper:
    def __init__(
        self,
        surrogate: SecondaryLossSurrogate | None = None,
        resistivity: ResistivityTable | None = None,
        pulse: TrapezoidPulse = TrapezoidPulse(),
        n_secondary=8,
        n_strands_total=2432,
        strand_diameter_m=0.2e-3,
        aluminum_density_kg_m3=2700.0,
    ):
        self.surrogate = SecondaryLossSurrogate() if surrogate is None else surrogate
        self.rho = ResistivityTable() if resistivity is None else resistivity
        self.pulse = pulse
        self.n_secondary = int(n_secondary)
        self.n_strands_total = int(n_strands_total)
        self.strand_diameter_m = float(strand_diameter_m)
        self.A_al_m2 = self.n_strands_total * np.pi * (self.strand_diameter_m / 2) ** 2
        self.mass_per_m = aluminum_density_kg_m3 * self.A_al_m2
        self.enthalpy = EnthalpyModel()
        self._build_C_splines()

    def _build_C_splines(self):
        db = self.surrogate.df
        self.T_nodes = self.surrogate.T_grid
        self.B_nodes = self.surrogate.B_grid
        self.C = {}
        for T in self.T_nodes:
            d = db[np.isclose(db["T_K"], T)].sort_values("Bpk_T")
            b = d["Bpk_T"].to_numpy(float)
            for col in ["Qxx_Jpm", "Qyy_Jpm", "Qxy_Jpm"]:
                q = d[col].to_numpy(float)
                F = np.zeros_like(b)
                F[1:] = q[1:] * self.pulse.rise_s / (2 * b[1:])
                self.C[(float(T), col)] = PchipInterpolator(b, F, extrapolate=True).derivative()

    def _C_component(self, B, T, col):
        B, T = np.broadcast_arrays(np.asarray(B, float), np.asarray(T, float))
        bf = np.clip(B.ravel(), self.B_nodes[0], self.B_nodes[-1])
        vals = np.stack([self.C[(float(Tn), col)](bf) for Tn in self.T_nodes], axis=0)
        out = np.empty(B.size)
        Tf = T.ravel()
        for i in range(B.size):
            out[i] = np.interp(Tf[i], self.T_nodes, vals[:, i])
        return out.reshape(B.shape)

    def run(
        self,
        field_model: FiniteTurnBiotSavart,
        target_order=2,
        dt_s=20e-6,
        initial_temperature_K=4.2,
        current_peak_A=47.29e3,
    ):
        s = field_model.turn_gauss_points(target_order)
        Bpk = s["Bmag_pk_T"].to_numpy(float)
        theta = np.arctan2(s["Bz_pk_T"], s["Br_pk_T"])
        c = np.cos(theta)
        ss = np.sin(theta)
        weights = s["weight"].to_numpy(float)
        length = s["turn_length_m"].to_numpy(float)

        Tloc = np.full(len(s), initial_temperature_K)
        hloc = np.zeros(len(s))
        Qeddy = np.zeros(len(s))
        Qj = np.zeros(len(s))
        times = np.arange(0.0, self.pulse.duration_s + 0.5 * dt_s, dt_s)
        history = []

        for k in range(len(times) - 1):
            tm = 0.5 * (times[k] + times[k + 1])
            g = self.pulse.amplitude(tm)
            dgdt = self.pulse.derivative(tm)
            Bnow = Bpk * g
            dBdt = Bpk * dgdt
            Cxx = self._C_component(Bnow, Tloc, "Qxx_Jpm")
            Cyy = self._C_component(Bnow, Tloc, "Qyy_Jpm")
            Cxy = self._C_component(Bnow, Tloc, "Qxy_Jpm")
            Ctheta = Cxx * c**2 + Cyy * ss**2 + 2 * Cxy * ss * c
            Peddy_per_m = self.n_secondary * Ctheta * dBdt**2

            I = current_peak_A * g
            rho = self.rho(Tloc, Bnow)
            Pj_per_m = I**2 * rho / self.A_al_m2

            Qeddy += Peddy_per_m * dt_s
            Qj += Pj_per_m * dt_s
            hloc += (Peddy_per_m + Pj_per_m) * dt_s / self.mass_per_m
            Tloc = self.enthalpy.temperature_from_specific_enthalpy(hloc)

            Peddy_total = float(np.sum(Peddy_per_m * weights * length))
            Pj_total = float(np.sum(Pj_per_m * weights * length))
            wmass = weights * length
            history.append(
                {
                    "time_s": tm,
                    "P_eddy_W": Peddy_total,
                    "P_transport_W": Pj_total,
                    "T_mean_K": float(np.average(Tloc, weights=wmass)),
                    "T_max_K": float(np.max(Tloc)),
                }
            )

        s["Qeddy_local_Jpm"] = Qeddy
        s["Qtransport_local_Jpm"] = Qj
        s["Tend_K"] = Tloc

        turn_rows = []
        for tid, gdf in s.groupby("turn_id"):
            wt = gdf["weight"].to_numpy(float)
            L = float(gdf["turn_length_m"].iloc[0])
            qe = float(np.sum(wt * gdf["Qeddy_local_Jpm"]))
            qj = float(np.sum(wt * gdf["Qtransport_local_Jpm"]))
            Tavg = float(np.sum(wt * gdf["Tend_K"]))
            row = gdf.iloc[0]
            turn_rows.append(
                {
                    "turn_id": int(tid),
                    "sp_id": int(row.sp_id),
                    "dp_id": int(row.dp_id),
                    "radial_layer": int(row.radial_layer),
                    "r_m": float(row.r_center_m),
                    "z_m": float(row.z_center_m),
                    "Qeddy_J": qe * L,
                    "Qtransport_J": qj * L,
                    "Qtotal_J": (qe + qj) * L,
                    "Tavg_end_K": Tavg,
                    "Tmax_end_K": float(gdf["Tend_K"].max()),
                }
            )
        turns = pd.DataFrame(turn_rows)
        hist = pd.DataFrame(history)
        summary = MagnetLossSummary(
            Q_eddy_J=float(turns["Qeddy_J"].sum()),
            Q_transport_J=float(turns["Qtransport_J"].sum()),
            Q_total_J=float(turns["Qtotal_J"].sum()),
            T_mean_end_K=float(np.average(s["Tend_K"], weights=weights * length)),
            T_max_end_K=float(s["Tend_K"].max()),
        )
        return s, turns, hist, summary
