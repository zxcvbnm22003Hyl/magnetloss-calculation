from apexloss import MagnetGeometry, FiniteTurnBiotSavart, SecondaryLossSurrogate
from apexloss.whole_magnet import fixed_temperature_loss

field = FiniteTurnBiotSavart(MagnetGeometry(), current_A=47.29e3, source_order=20)
surr = SecondaryLossSurrogate()
turns, Q = fixed_temperature_loss(field, surr, temperature_K=4.2, target_order=2)
print("whole-magnet intrinsic strand eddy loss [kJ/pulse] =", Q / 1e3)
print(turns.sort_values("Q_eddy_J", ascending=False).head())
