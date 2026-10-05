from apexloss import MagnetGeometry, FiniteTurnBiotSavart
from apexloss.whole_magnet import ModalThermalMapper

field = FiniteTurnBiotSavart(
    MagnetGeometry(),
    current_A=47.29e3,
    source_order=20,
)

mapper = ModalThermalMapper(n_modes=12)

samples, turns, history, summary = mapper.run(
    field,
    target_order=2,
    dt_s=5e-6,
    initial_temperature_K=4.2,
    current_peak_A=47.29e3,
)

print(summary)
