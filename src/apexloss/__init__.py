"""APEX-PulseLoss: reduced-order electromagnetic loss models for pulsed magnets."""

from .geometry import MagnetGeometry
from .materials import ResistivityTable, EnthalpyModel
from .waveforms import TrapezoidPulse
from .field import FiniteTurnBiotSavart
from .strand_diffusion import SingleStrandDiffusion
from .modal_diffusion import RobinBesselStrand
from .surrogate import SecondaryLossSurrogate

__all__ = [
    "MagnetGeometry",
    "ResistivityTable",
    "EnthalpyModel",
    "TrapezoidPulse",
    "FiniteTurnBiotSavart",
    "SingleStrandDiffusion",
    "RobinBesselStrand",
    "SecondaryLossSurrogate",
]

__version__ = "0.1.1"
