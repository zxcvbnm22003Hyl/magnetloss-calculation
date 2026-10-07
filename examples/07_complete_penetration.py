from apexloss import (
    CompletePenetrationStrand,
    FiniteTurnBiotSavart,
    MagnetGeometry,
)
from apexloss.materials import LegacyKohlerRRR
from apexloss.whole_magnet import CompletePenetrationThermalMapper

# 1) Single-strand analytical baseline
strand = CompletePenetrationStrand(diameter_m=0.2e-3)
single = strand.simulate_trapezoid(
    B_peak_T=20.619192,
    dt_s=1e-6,
    resistivity=LegacyKohlerRRR(),
)
print("single-strand CP:", single)

# 2) Whole-magnet rho(B,T) + adiabatic feedback
field = FiniteTurnBiotSavart(
    MagnetGeometry(),
    current_A=47.29e3,
    source_order=20,
)

mapper = CompletePenetrationThermalMapper()
samples, turns, history, summary = mapper.run(
    field,
    target_order=2,
    dt_s=5e-6,
    initial_temperature_K=4.2,
    current_peak_A=47.29e3,
)
print("whole-magnet CP:", summary)
