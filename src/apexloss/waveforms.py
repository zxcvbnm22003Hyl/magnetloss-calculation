from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class TrapezoidPulse:
    """Normalized trapezoidal waveform with linear rise, flat top and linear decay."""

    rise_s: float = 3e-3
    flat_s: float = 1e-3
    decay_s: float = 3e-3

    @property
    def duration_s(self) -> float:
        return self.rise_s + self.flat_s + self.decay_s

    def amplitude(self, t):
        t = np.asarray(t, dtype=float)
        y = np.zeros_like(t)
        m = (t >= 0.0) & (t < self.rise_s)
        y[m] = t[m] / self.rise_s
        m = (t >= self.rise_s) & (t < self.rise_s + self.flat_s)
        y[m] = 1.0
        m = (t >= self.rise_s + self.flat_s) & (t <= self.duration_s)
        y[m] = (self.duration_s - t[m]) / self.decay_s
        return float(y) if y.ndim == 0 else y

    def derivative(self, t):
        t = np.asarray(t, dtype=float)
        y = np.zeros_like(t)
        y[(t >= 0.0) & (t < self.rise_s)] = 1.0 / self.rise_s
        y[(t >= self.rise_s + self.flat_s) & (t <= self.duration_s)] = -1.0 / self.decay_s
        return float(y) if y.ndim == 0 else y
