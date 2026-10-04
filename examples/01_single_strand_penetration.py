from apexloss import SingleStrandDiffusion
from apexloss.materials import LegacyKohlerRRR

model = SingleStrandDiffusion(diameter_m=0.2e-3, radial_cells=160)
result = model.simulate_trapezoid(
    B_peak_T=20.619192,
    dt_s=1e-6,
    tail_s=23e-3,
    resistivity=LegacyKohlerRRR(),
)

print("Finite-diffusion pulse loss [J/m] :", result.Q_pulse_J_per_m)
print("Post-pulse tail [J/m]             :", result.Q_tail_J_per_m)
print("Total finite-diffusion [J/m]      :", result.Q_total_J_per_m)
print("Complete-penetration [J/m]        :", result.Q_complete_penetration_J_per_m)
