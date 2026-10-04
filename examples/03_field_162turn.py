from apexloss import MagnetGeometry, FiniteTurnBiotSavart

field = FiniteTurnBiotSavart(MagnetGeometry(), current_A=47.29e3, source_order=20)
print("B0 [T] =", field.center_field_T())
turns = field.turn_center_fields()
print(turns.head())
print("max turn-center |B| [T] =", turns.Bmag_pk_T.max())
