# Changelog

## 0.1.2 — 2026-10-07

- add dedicated complete-penetration analytical strand solver
- add `whole-magnet --mode cp-thermal` with rho(B,T) and adiabatic feedback
- add `strand-cp` CLI command
- add CP regression tests, example and frozen reference results
- centralize the Robin–Bessel fast-diffusion limit on the CP implementation


## 0.1.1 — 2026-10-05

- add Robin–Bessel finite-penetration modal strand solver
- add state-resolved whole-magnet modal thermal mapper
- expose `whole-magnet --mode modal-thermal`
- validate modal strand loss against the radial diffusion reference
- add hotspot 2432-strand FEM and whole-magnet comparison references


## 0.1.0 — 2026-10-04

Initial public-code package structure.

- finite rectangular-turn Biot-Savart solver
- single-strand finite magnetic penetration / diffusion solver
- explicit 1/19/304-strand 2-D P1 FEM
- bundled 304-strand loss surrogate
- 162-turn fixed-temperature and adiabatic reduced-order loss mapping
- legacy pure-Python scripts retained for provenance
- regression tests and reference outputs
