from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .materials import ResistivityTable
from .waveforms import TrapezoidPulse


@dataclass
class CompletePenetrationResult:
    """Single-strand complete-penetration loss result."""

    Q_pulse_J_per_m: float
    Q_tail_J_per_m: float
    Q_total_J_per_m: float
    P_peak_W_per_m: float
    t_P_peak_s: float


class CompletePenetrationStrand:
    """Complete-penetration eddy-loss model for a round normal-metal strand.

    For a round strand of radius a in a spatially uniform transverse magnetic
    field B(t), the complete-penetration approximation gives

        P'/L = pi a^4 / (4 rho) * (dB/dt)^2.

    The model is memoryless: the instantaneous loss depends only on the current
    rho(B,T) and dB/dt. Therefore the eddy-current loss is zero whenever
    dB/dt = 0, including the flat top and any post-pulse tail.

    This class is retained as the analytical baseline for comparison with the
    finite-diffusion Robin-Bessel and radial-diffusion models.
    """

    def __init__(self, diameter_m=0.2e-3):
        self.diameter_m = float(diameter_m)
        self.radius_m = 0.5 * self.diameter_m

    @staticmethod
    def power_W_per_m(radius_m, rho_ohm_m, dBdt_T_s):
        """Instantaneous complete-penetration Joule power per unit length."""
        radius = float(radius_m)
        rho = np.asarray(rho_ohm_m, dtype=float)
        dBdt = np.asarray(dBdt_T_s, dtype=float)
        return np.pi * radius**4 / (4.0 * rho) * dBdt**2

    def power(self, rho_ohm_m, dBdt_T_s):
        """Instance form of :meth:`power_W_per_m`."""
        return self.power_W_per_m(self.radius_m, rho_ohm_m, dBdt_T_s)

    def simulate_trapezoid(
        self,
        B_peak_T=20.619192,
        pulse: TrapezoidPulse = TrapezoidPulse(),
        dt_s=1e-6,
        tail_s=0.0,
        temperature_K=4.2,
        resistivity=None,
    ) -> CompletePenetrationResult:
        """Integrate the complete-penetration loss over a trapezoidal pulse.

        Parameters
        ----------
        B_peak_T:
            Peak transverse magnetic field.
        pulse:
            Normalized drive waveform.
        dt_s:
            Midpoint-integration time step.
        tail_s:
            Accepted for API symmetry with finite-diffusion models. The CP
            model has no magnetic state memory, so the post-pulse eddy loss is
            identically zero.
        temperature_K:
            Fixed strand temperature for this single-strand calculation.
        resistivity:
            Either ResistivityTable(T,B), a callable rho(B), or None to use the
            bundled ResistivityTable.
        """
        table = ResistivityTable() if resistivity is None else resistivity

        def rho_of_B(B):
            if isinstance(table, ResistivityTable):
                return float(table(temperature_K, B))
            return float(table(B))

        Qpulse = 0.0
        Pmax = 0.0
        tPmax = 0.0

        nsteps = int(round(pulse.duration_s / dt_s))
        for k in range(nsteps):
            tm = (k + 0.5) * dt_s
            B = B_peak_T * pulse.amplitude(tm)
            dBdt = B_peak_T * pulse.derivative(tm)
            rho = rho_of_B(B)
            power = float(self.power(rho, dBdt))
            Qpulse += power * dt_s
            if power > Pmax:
                Pmax = power
                tPmax = tm

        return CompletePenetrationResult(
            Q_pulse_J_per_m=Qpulse,
            Q_tail_J_per_m=0.0,
            Q_total_J_per_m=Qpulse,
            P_peak_W_per_m=Pmax,
            t_P_peak_s=tPmax,
        )

    def closed_form_constant_rho(
        self,
        B_peak_T,
        rho_ohm_m,
        pulse: TrapezoidPulse = TrapezoidPulse(),
    ):
        """Closed-form CP energy for constant rho and a trapezoidal pulse."""
        return (
            np.pi * self.radius_m**4
            / (4.0 * float(rho_ohm_m))
            * float(B_peak_T) ** 2
            * (1.0 / pulse.rise_s + 1.0 / pulse.decay_s)
        )
