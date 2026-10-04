# Numerical methods

## Finite rectangular-turn field integration

A turn occupying a rectangular region in `(r,z)` is decomposed into differential circular loops. Gauss-Legendre quadrature integrates the circular-loop field kernel over the source cross section.

An even source quadrature order is used to avoid placing a source quadrature node exactly at a turn center during self-field evaluation.

## Single-strand finite penetration

The radial `m=1` magnetic-diffusion operator is discretized using the original APEX finite-difference formulation. Crank-Nicolson time stepping leads to a tridiagonal linear system solved with `scipy.linalg.solve_banded`.

## 2-D multifilament FEM

The scalar magnetic vector potential `A_z` is discretized with linear triangular basis functions. The conductor mass matrix is projected to remove one scalar mean electric-field mode per strand. This enforces zero net current independently in every insulated strand while retaining local eddy-current loops.

## Loss surrogate

For fixed 3/1/3 ms waveform, the 304-strand full-pulse response is tabulated versus peak field, field direction, and temperature. Low-field interpolation is performed on `Q/B^2`, not directly on `Q`, to preserve the correct quadratic field scaling.

## Fast thermal reconstruction

The reduced-order thermal model assumes

```math
P_{eddy}=C(B,T,\theta)\left(\frac{dB}{dt}\right)^2.
```

For a symmetric linear ramp, `C` is recovered by differentiating the integral relation implied by the full-pulse surrogate. This is computationally inexpensive but does not preserve the full internal magnetic-diffusion state.
