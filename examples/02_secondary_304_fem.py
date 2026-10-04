from apexloss.multifilament_fem import MultifilamentFEM

fem = MultifilamentFEM.secondary304(mesh_h_m=0.03e-3, domain_half_m=2.8e-3)
result = fem.run_trapezoid(
    B_peak_T=20.619192,
    angle_rad=1.5707963267948966,
    temperature_K=4.2,
    dt_s=0.25e-3,
)

print("304-strand field-only loss [J/m] :", result.energy_J_per_m)
print("nodes                            :", result.n_nodes)
print("triangles                        :", result.n_triangles)
print("mean strand area error          :", result.mean_strand_area_error)
