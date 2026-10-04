from apexloss import MagnetGeometry, FiniteTurnBiotSavart
from apexloss.whole_magnet import QuasiStaticThermalMapper

field = FiniteTurnBiotSavart(MagnetGeometry(), current_A=47.29e3, source_order=20)
mapper = QuasiStaticThermalMapper()
samples, turns, history, summary = mapper.run(
    field, target_order=2, dt_s=20e-6, initial_temperature_K=4.2
)
print(summary)
