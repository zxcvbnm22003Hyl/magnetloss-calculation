# Validation snapshot for v0.1.0

The following values are frozen as regression references. They are not uncertainty-free engineering truth; they are numerical consistency targets for this release.

## Single 0.2 mm strand finite penetration

Conditions:

- Bpk = 20.619192 T
- 3/1/3 ms trapezoid
- 4.2 K legacy Kohler/RRR resistivity
- 160 radial cells
- dt = 1 us
- diffusion tail integrated to 30 ms total time

Results:

```text
Q(0-7 ms)                 = 0.5810907987 J/m
Q(7-30 ms tail)           = 0.0072395840 J/m
Q(total finite diffusion) = 0.5883303827 J/m
Q(complete penetration)   = 0.6083368195 J/m
```

Thus the finite-diffusion result is below the complete-penetration reference, as expected.

## 162-turn finite-cross-section field

At Ipk = 47.29 kA and source quadrature order 20:

```text
B0 = 20.195228537 T
```

The maximum turn-center field occurs at the innermost turns adjacent to the mid-plane and is approximately 20.63 T.

## 304-strand surrogate, fixed 4.2 K

Using per-turn 2x2 Gauss mapping, the current v0.1.0 Python package gives approximately:

```text
Qeddy ≈ 73.34 kJ/pulse
```

The 2x2, 3x3 and 4x4 macro cross-section integrations are effectively converged relative to the local-model uncertainty.

## Thermal reduced-order mapping

The fast surrogate-based adiabatic model gives a whole-magnet intrinsic eddy loss around 19.9 kJ/pulse and transport Joule heat around 2.3 kJ/pulse for the current baseline.

An independent older 648-point single-strand spatial-RVE model gives 18.8466 kJ/pulse eddy loss. The difference is retained as model-form disagreement rather than removed by calibration.


## Robin-Bessel modal validation

Using the same 0.2 mm strand, 20.619192 T peak field and 3/1/3 ms trapezoid as the radial diffusion regression:

```text
20-mode Robin-Bessel, dt=1 us:
Q(0-7 ms)   = 0.5810597 J/m
Q(7-30 ms)  = 0.00724007 J/m
Q(total)    = 0.5882998 J/m

160-cell radial diffusion:
Q(total)    = 0.5883304 J/m
```

The total-energy difference is approximately -0.005%.

For the dynamic hotspot tertiary cable with the same rho(B,T) and adiabatic enthalpy model:

```text
complete penetration     = 276.139 J/m
Robin-Bessel modal       = 274.718 J/m
explicit 2432-strand FEM = 273.749 J/m
```

The modal result is about 0.35% above the explicit 2432-strand FEM result.

## Whole-magnet dynamic comparison

For the 162-turn APEX reference magnet, 2x2 Gauss sampling per turn and the RRR3000 rho(B,T) model:

```text
complete penetration + thermal feedback = 20.959 kJ/pulse
Robin-Bessel modal + thermal feedback    = 20.707 kJ/pulse
old 648-point spatial-RVE                = 21.064 kJ/pulse
304/2432 reduced-order surrogate chain   = 22.958 kJ/pulse
```

The first three independent state-resolved/strand-resolved approaches cluster around 20.7--21.1 kJ/pulse. The higher surrogate-chain value is retained as model-form disagreement and should not be removed by calibration; its dynamic state reconstruction remains under review.


## Field-strength validation scan

The dynamic rho(B,T) model was checked at 5, 10, 15 and 20.5 T. The explicit
304-strand FEM used 31.25 and 62.5 us time steps and was first-order
Richardson-extrapolated to dt -> 0. The Robin-Bessel model used 12 modes and a
2 us time step.

| Bpk (T) | FEM304 total (J/m) | modal total (J/m) | modal - FEM |
|---:|---:|---:|---:|
| 5.0  | 7.446686  | 7.405439  | -0.554% |
| 10.0 | 16.488514 | 16.455709 | -0.199% |
| 15.0 | 25.057773 | 25.049136 | -0.034% |
| 20.5 | 34.134005 | 34.150586 | +0.049% |

Thus the maximum absolute modal/FEM discrepancy across 5--20.5 T is about
0.56%, with no systematic high-field divergence.

A direct 2432-strand / 304-strand paired calculation gives

| Bpk (T) | C3_total = Q2432 / (8 Q304) |
|---:|---:|
| 5.0  | 1.001405 |
| 10.0 | 1.000381 |
| 15.0 | 1.000270 |
| 20.5 | 1.000628 |

With matched conductor quadrature and discretization, the tertiary magnetic
interaction correction is therefore only about 0.03--0.14% over the scanned
field range. The larger tertiary corrections obtained in earlier exploratory
runs were dominated by mismatched conductor-integration / mesh formulations.

Applying C3(B) to the 162-turn Robin-Bessel whole-magnet result changes the
total pulse loss only from about 20.717 kJ to 20.729 kJ, i.e. +0.055%.
The recommended current whole-magnet value is therefore

```text
Qmag ≈ 20.73 kJ/pulse
```

under the present insulated-strand, adiabatic, rho(B,T) model assumptions.
