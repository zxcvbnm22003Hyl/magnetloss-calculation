# Model assumptions

## Electromagnetic decomposition

The conductor is assumed to contain electrically insulated strands. The default loss budget is

```math
Q_{\rm total}\approx Q_{\rm transport}+Q_{\rm intra-filament}.
```

Inter-strand coupling loss is excluded:

```math
Q_{\rm coupling}\approx 0.
```

This assumption must be revisited when measured transverse/contact resistance becomes available.

## Single-strand diffusion model

The strand is circular and long relative to its diameter. The external field is transverse to the strand axis. The numerical model solves the `m=1` transverse magnetic-diffusion mode in radius and retains finite magnetic penetration.

The complete-penetration formula is provided only as a reference model.

## Multifilament FEM

The explicit 2-D model assumes straight, untwisted strands. Each strand is electrically isolated and has zero net transport current in field-only calculations. The current repository does not contain a helicoidal coordinate transformation.

The geometry is represented on a structured triangular background mesh. Strand boundaries are therefore stair-stepped/embedded rather than body-fitted.

## Whole-magnet field

The magnet is air-cored and axisymmetric. Each rectangular macro turn carries uniform azimuthal current density and is integrated as circular loop quadrature.

Because the magnetic problem is linear in current, a common-drive waveform gives

```math
B_r(t)=g(t)B_{r,pk},\qquad B_z(t)=g(t)B_{z,pk}.
```

## Whole-magnet loss mapping

The default whole-magnet mapping uses 2x2 Gauss points per macro-turn cross section. For the current geometry, 2x2, 3x3 and 4x4 mappings are numerically converged in total loss to substantially below the local-model uncertainty.

## Thermal model

The thermal model is adiabatic for one pulse. It uses the existing aluminum heat-capacity proxy and does not model coolant, insulation, metallic reinforcement, contact resistance, or heat transfer between turns.

## Material data

The bundled resistivity table is a design-level `rho(T,B)` model for high-purity aluminum with an effective RRR of 3000. It must eventually be replaced by measured finished-strand/cable data.
