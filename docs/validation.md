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
