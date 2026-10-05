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


## Robin-Bessel finite-penetration modal model

For a round strand in a transverse field, write the internal vector potential as

```math
A_z(r,\phi,t)=f(r,t)\sin\phi.
```

Matching the interior solution to the exterior vacuum field including the induced dipole term gives the Robin boundary condition

```math
f(a)+af'(a)=2aB_{\rm ext}.
```

After subtracting the imposed field contribution, the homogeneous radial eigenfunctions are `J1(lambda_n r/a)`, where

```math
J_0(\lambda_n)=0.
```

The modal amplitudes satisfy

```math
\dot q_n+q_n/\tau_n=-\beta_n\dot B,
```

with

```math
\tau_n=\mu_0a^2/[\rho(B,T)\lambda_n^2],
\qquad
\beta_n=4/[\lambda_n^2J_1(\lambda_n)].
```

The instantaneous strand loss is

```math
P'=
\frac{\pi a^4}{2\rho}
\sum_nJ_1^2(\lambda_n)(q_n/\tau_n)^2.
```

This formulation retains magnetic-diffusion memory through the rise, flat top, decay and optional post-pulse tail. The modal ODE is integrated analytically over each time step with midpoint-frozen `rho(B,T)` and `dB/dt`. In the fast-diffusion limit it recovers the complete-penetration formula exactly.

For the present 0.2 mm strand, 10--12 modes are already sufficient for engineering calculations; 20 modes are used for the frozen single-strand regression.
