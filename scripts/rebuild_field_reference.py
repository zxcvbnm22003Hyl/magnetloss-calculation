from apexloss import MagnetGeometry, FiniteTurnBiotSavart

field = FiniteTurnBiotSavart(MagnetGeometry(), current_A=47.29e3, source_order=20)
turns = field.turn_center_fields()
turns.to_csv("APEX_162turn_peak_fields_finite_rect.csv", index=False)
print("B0 [T] =", field.center_field_T())
print("Bmax turn center [T] =", turns.Bmag_pk_T.max())
