from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.special import jn_zeros, jv

from .constants import MU0
from .materials import ResistivityTable
from .waveforms import TrapezoidPulse


@dataclass
class ModalStrandResult:
    Q_pulse_J_per_m: float
    Q_tail_J_per_m: float
    Q_total_J_per_m: float
    P_peak_W_per_m: float
    t_P_peak_s: float


class RobinBesselStrand:
    """Finite-penetration modal model for a round strand in transverse field.

    The external vacuum field is matched with the Robin boundary condition

        f(a) + a f'(a) = 2 a B_ext,

    so the radial eigenvalues are the positive zeros of J0, not J1.

    The modal state q_n obeys

        dq_n/dt + q_n/tau_n = -beta_n dB/dt,

    with

        tau_n  = mu0 a^2 / (rho lambda_n^2)
        beta_n = 4 / (lambda_n^2 J1(lambda_n)).

    The instantaneous Joule loss per unit length is

        P' = pi a^4/(2 rho) *
             sum_n J1(lambda_n)^2 (q_n/tau_n)^2.

    In the fast-diffusion limit this reduces exactly to the
    complete-penetration formula pi a^4/(4 rho) (dB/dt)^2.
    """

    def __init__(self, diameter_m=0.2e-3, n_modes=12):
        self.diameter_m = float(diameter_m)
        self.radius_m = 0.5 * self.diameter_m
        self.n_modes = int(n_modes)
        if self.n_modes < 1:
            raise ValueError("n_modes must be >= 1")
        self.lambda_n = jn_zeros(0, self.n_modes)
        self.J1_n = jv(1, self.lambda_n)
        self.beta_n = 4.0 / (self.lambda_n**2 * self.J1_n)

    def initial_state(self, shape=()):
        """Return zero modal state with trailing modal dimension."""
        return np.zeros(tuple(shape) + (self.n_modes,), dtype=float)

    def time_constants_s(self, rho_ohm_m):
        rho = np.asarray(rho_ohm_m, dtype=float)
        return (
            MU0 * self.radius_m**2
            / (rho[..., None] * self.lambda_n**2)
        )

    def step(self, state, dt_s, dBdt_T_s, rho_ohm_m):
        """Advance the modal state over one time step.

        rho and dB/dt are frozen at the time-step midpoint. The linear modal
        ODE is then integrated analytically over the step.

        Parameters
        ----------
        state:
            Array with trailing dimension n_modes.
        dt_s:
            Time-step size.
        dBdt_T_s:
            Scalar or array broadcastable to state.shape[:-1].
        rho_ohm_m:
            Scalar or array broadcastable to state.shape[:-1].

        Returns
        -------
        state_new, power_W_per_m
        """
        q = np.asarray(state, dtype=float)
        if q.shape[-1] != self.n_modes:
            raise ValueError("state trailing dimension must equal n_modes")

        rho = np.asarray(rho_ohm_m, dtype=float)
        dBdt = np.asarray(dBdt_T_s, dtype=float)
        tau = self.time_constants_s(rho)

        e = np.exp(-dt_s / tau)
        eh = np.exp(-0.5 * dt_s / tau)
        forcing = self.beta_n * dBdt[..., None] * tau

        q_mid = q * eh - forcing * (1.0 - eh)
        q_new = q * e - forcing * (1.0 - e)

        power = (
            np.pi * self.radius_m**4 / (2.0 * rho)
            * np.sum(
                self.J1_n**2 * (q_mid / tau) ** 2,
                axis=-1,
            )
        )
        return q_new, power

    @staticmethod
    def complete_penetration_power_W_per_m(radius_m, rho_ohm_m, dBdt_T_s):
        return (
            np.pi * float(radius_m) ** 4
            / (4.0 * np.asarray(rho_ohm_m, dtype=float))
            * np.asarray(dBdt_T_s, dtype=float) ** 2
        )

    def simulate_trapezoid(
        self,
        B_peak_T=20.619192,
        pulse: TrapezoidPulse = TrapezoidPulse(),
        dt_s=1e-6,
        tail_s=23e-3,
        temperature_K=4.2,
        resistivity=None,
    ) -> ModalStrandResult:
        """Fixed-temperature single-strand regression/validation calculation."""
        table = ResistivityTable() if resistivity is None else resistivity

        def rho_of_B(B):
            if isinstance(table, ResistivityTable):
                return float(table(temperature_K, B))
            return float(table(B))

        q = self.initial_state()
        Qpulse = 0.0
        Qtail = 0.0
        Pmax = 0.0
        tPmax = 0.0
        tend = pulse.duration_s + tail_s
        nsteps = int(round(tend / dt_s))

        for k in range(nsteps):
            tm = (k + 0.5) * dt_s
            amp = pulse.amplitude(tm)
            B = B_peak_T * amp
            dBdt = B_peak_T * pulse.derivative(tm)
            rho = rho_of_B(B)
            q, power = self.step(q, dt_s, dBdt, rho)
            power = float(power)
            if tm <= pulse.duration_s:
                Qpulse += power * dt_s
            else:
                Qtail += power * dt_s
            if power > Pmax:
                Pmax = power
                tPmax = tm

        return ModalStrandResult(
            Q_pulse_J_per_m=Qpulse,
            Q_tail_J_per_m=Qtail,
            Q_total_J_per_m=Qpulse + Qtail,
            P_peak_W_per_m=Pmax,
            t_P_peak_s=tPmax,
        )
