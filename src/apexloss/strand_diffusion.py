from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import solve_banded

from .constants import MU0
from .materials import ResistivityTable, LegacyKohlerRRR
from .waveforms import TrapezoidPulse


@dataclass
class StrandDiffusionResult:
    Q_pulse_J_per_m: float
    Q_tail_J_per_m: float
    Q_total_J_per_m: float
    P_peak_W_per_m: float
    t_P_peak_s: float
    Q_complete_penetration_J_per_m: float | None = None


class SingleStrandDiffusion:
    """Finite magnetic penetration model for one round strand in transverse field."""

    def __init__(self, diameter_m=0.2e-3, radial_cells=160):
        self.diameter_m = float(diameter_m)
        self.radius_m = self.diameter_m / 2
        self.radial_cells = int(radial_cells)
        self._build_operator()

    def _build_operator(self):
        N = self.radial_cells
        a = self.radius_m
        dr = a / N
        r = np.arange(N) * dr
        lo = np.zeros(N)
        dg = np.zeros(N)
        up = np.zeros(N)
        b = np.zeros(N)
        dg[0] = -8 / dr**2
        up[0] = 8 / dr**2
        q = a / (2 * dr)
        den = 2 + 3 * q
        alphaB = 2 / den
        c1 = 4 * q / den
        c2 = -q / den
        for i in range(1, N):
            ri = r[i]
            cm = 1 / dr**2 - 3 / (2 * ri * dr)
            c0 = -2 / dr**2
            cp = 1 / dr**2 + 3 / (2 * ri * dr)
            lo[i] += cm
            dg[i] += c0
            if i < N - 1:
                up[i] += cp
            else:
                dg[i] += cp * c1
                lo[i] += cp * c2
                b[i] += cp * alphaB
        self.r, self.dr, self.lo, self.dg, self.up, self.bc = r, dr, lo, dg, up, b

    def _tri_mv(self, x):
        y = self.dg * x
        y[1:] += self.lo[1:] * x[:-1]
        y[:-1] += self.up[:-1] * x[1:]
        return y

    def simulate_trapezoid(
        self,
        B_peak_T=20.619192,
        pulse: TrapezoidPulse = TrapezoidPulse(),
        dt_s=1e-6,
        tail_s=23e-3,
        temperature_K=4.2,
        resistivity=None,
        complete_penetration_reference=True,
    ) -> StrandDiffusionResult:
        if resistivity is None:
            table = ResistivityTable()
            rho_fun = lambda B: float(table(temperature_K, B))
        elif isinstance(resistivity, ResistivityTable):
            rho_fun = lambda B: float(resistivity(temperature_K, B))
        else:
            rho_fun = resistivity

        a = self.radius_m
        N = self.radial_cells
        tend = pulse.duration_s + tail_s
        nsteps = int(round(tend / dt_s))
        g = np.zeros(N)
        Qpulse = 0.0
        Qtail = 0.0
        Pmax = 0.0
        tPmax = 0.0

        for n in range(nsteps):
            tm = (n + 0.5) * dt_s
            Bm = B_peak_T * pulse.amplitude(tm)
            rho = rho_fun(Bm)
            D = rho / MU0
            c = 0.5 * dt_s * D
            rhs = g + c * self._tri_mv(g) + dt_s * D * self.bc * Bm
            ab = np.zeros((3, N))
            ab[0, 1:] = -c * self.up[:-1]
            ab[1, :] = 1 - c * self.dg
            ab[2, :-1] = -c * self.lo[1:]
            gnew = solve_banded((1, 1), ab, rhs)
            gmid = 0.5 * (g + gnew)
            Lg = self._tri_mv(gmid) + self.bc * Bm
            Jamp = -(self.r / MU0) * Lg
            Jb = 2 * Jamp[-1] - Jamp[-2]
            rr = np.r_[self.r, a]
            JJ = np.r_[Jamp, Jb]
            P = np.pi * rho * np.trapezoid(JJ**2 * rr, rr)
            if tm <= pulse.duration_s:
                Qpulse += P * dt_s
            else:
                Qtail += P * dt_s
            if P > Pmax:
                Pmax = P
                tPmax = tm
            g = gnew

        Qcp = None
        if complete_penetration_reference:
            Qcp = self.complete_penetration_loss(B_peak_T, pulse, temperature_K, resistivity)
        return StrandDiffusionResult(Qpulse, Qtail, Qpulse + Qtail, Pmax, tPmax, Qcp)

    def complete_penetration_loss(
        self,
        B_peak_T,
        pulse=TrapezoidPulse(),
        temperature_K=4.2,
        resistivity=None,
        dt_s=1e-7,
    ):
        if resistivity is None:
            table = ResistivityTable()
            rho_fun = lambda B: float(table(temperature_K, B))
        elif isinstance(resistivity, ResistivityTable):
            rho_fun = lambda B: float(resistivity(temperature_K, B))
        else:
            rho_fun = resistivity
        tt = np.arange(dt_s / 2, pulse.duration_s, dt_s)
        amp = pulse.amplitude(tt)
        B = B_peak_T * amp
        dBdt = B_peak_T * pulse.derivative(tt)
        rho = np.array([rho_fun(x) for x in B])
        return float(np.sum((np.pi * self.radius_m**4 / (4 * rho)) * dBdt**2 * dt_s))

    @staticmethod
    def legacy_default():
        return SingleStrandDiffusion(diameter_m=0.2e-3, radial_cells=160), LegacyKohlerRRR()
