# magnetloss-calculation

Pure-Python tools for magnetic-field calculation and pulsed-loss estimation in high-field solenoids using electrically insulated multifilament metallic conductors.

**Author:** hery rainforeset

The repository was developed around the APEX Phase-I pulsed-solenoid reference case, but the numerical components are separated so that the field solver, single-strand diffusion model, explicit multifilament FEM, surrogate model, and whole-magnet mapping can be used independently.

## What is included

The code contains four model levels:

1. **Finite-cross-section Biot–Savart solver** for axisymmetric solenoids.
2. **Single-strand finite magnetic-penetration model** for a round metallic strand in transverse pulsed field.
3. **Explicit 2-D P1 multifilament FEM** for 1, 19, or 304 electrically insulated strands.
4. **304-strand loss surrogate + whole-magnet mapping** with Gauss sampling and optional adiabatic thermal feedback.

The current package is entirely Python-based and does not require proprietary multiphysics software.

## Reference magnet

The default reference geometry is:

| Parameter | Value |
|---|---:|
| Double pancakes | 9 |
| Single pancakes | 18 |
| Radial turns per SP | 9 |
| Total turns | 162 |
| Inner radius | 125 mm |
| Outer radius | 215.7 mm |
| Active axial length | 332.7 mm |
| Cable envelope | 8.3 mm × 16.9 mm |
| Peak current | 47.29 kA |
| Pulse | 3 ms rise / 1 ms flat / 3 ms decay |
| Strand diameter | 0.2 mm |
| Strands per turn | 2432 |
| Effective RRR model | 3000 |

The modeled conductor hierarchy is

```text
0.2 mm strand
    -> 19-strand primary subcable
    -> 16 primary subcables = 304-strand secondary subcable
    -> 8 secondary subcables = 2432-strand turn cable
```

The default model assumes electrical insulation between strands and therefore excludes inter-strand coupling-current loss.

## Model hierarchy

### 1. Finite rectangular-turn Biot–Savart field solver

Each macroscopic turn is represented as a rectangular current-density region. Its field is integrated as a continuum of circular current loops using tensor-product Gauss–Legendre quadrature.

This avoids the self-field singularity that appears if a turn is represented by a zero-thickness filament and the field is evaluated at the turn center.

Current reference result at 47.29 kA:

```text
Center field B0                  ≈ 20.19523 T
Maximum turn-center |B|          ≈ 20.63009 T
```

Main implementation:

```text
src/apexloss/field.py
```

### 2. Single-strand finite magnetic penetration

`SingleStrandDiffusion` solves the time-domain magnetic-diffusion problem for the transverse `m = 1` mode of a round strand.

This model retains finite penetration and the post-pulse diffusion tail. It is independent of the 304-strand FEM and is intended as a numerical benchmark against the complete-penetration approximation

```math
Q_{\mathrm{CP}}
=
\int
\frac{\pi a^4}{4\rho(B,T)}
\left(\frac{dB}{dt}\right)^2
dt .
```

Frozen regression case:

```text
d                         = 0.2 mm
Bpk                       = 20.619192 T
pulse                     = 3/1/3 ms
radial cells              = 160
time step                 = 1 us

Q(0-7 ms)                 = 0.5810908 J/m
Q(post-pulse tail)        = 0.00723958 J/m
Q(total finite diffusion) = 0.5883304 J/m
Q(complete penetration)   = 0.6083368 J/m
```

Main implementation:

```text
src/apexloss/strand_diffusion.py
```

The original standalone numerical script is retained in

```text
legacy/APEX_strand_dynamic_rho_diffusion.py
```

for provenance and result traceability.

### 3. Explicit 2-D multifilament FEM

The local multifilament model uses a scalar magnetic vector potential `A_z` with linear triangular P1 elements.

Available conductor layouts:

- one 0.2 mm strand;
- one 19-strand `1+6+12` primary subcable;
- one 304-strand secondary subcable.

Each strand is electrically isolated. In field-only calculations, the strand-average electric-field mode is projected out independently for every strand, enforcing zero net transport current while retaining strand-internal eddy-current loops.

Main implementation:

```text
src/apexloss/multifilament_fem.py
```

The current mesh is an embedded structured triangular mesh rather than a body-fitted mesh. The remaining geometric discretization error should therefore be treated as part of the model uncertainty.

### 4. 304-strand surrogate and whole-magnet mapping

The bundled surrogate stores the full-pulse 304-strand response versus peak field and temperature.

Field-direction dependence is represented as

```math
Q(\theta)
=
Q_{xx}\cos^2\theta
+
Q_{yy}\sin^2\theta
+
2Q_{xy}\sin\theta\cos\theta .
```

For the full magnet, each turn cross section is sampled using Gauss points. The default is `2 × 2`, giving

```text
162 turns × 4 Gauss points = 648 macro field points
```

The whole-magnet calculation is therefore

```text
finite-turn Biot–Savart
    -> local Br, Bz at 648 Gauss points
    -> 304-strand surrogate
    -> ×8 secondary subcables
    -> integrate over each turn
    -> total magnet loss
```

For the current geometry, 2×2, 3×3, and 4×4 turn-cross-section sampling are already converged relative to the uncertainty of the local loss model.

## Current reference results

These values are regression targets for the current release, not uncertainty-free engineering truth.

| Quantity | Current reference |
|---|---:|
| Center field | 20.19523 T |
| Max turn-center field | 20.63009 T |
| Whole-magnet fixed-4.2 K intrinsic strand eddy loss | ~73.4 kJ/pulse |
| Older 648-point single-strand spatial-RVE eddy loss with thermal feedback | 18.8466 kJ/pulse |
| 304-strand reduced-order thermal result | ~19.9–20.3 kJ/pulse |
| Transport Joule heat in the current thermal model | ~2.3 kJ/pulse |
| End mean temperature | ~35.8 K |
| End local maximum temperature | ~47 K |

The disagreement between independent local models is kept explicitly as model-form uncertainty rather than removed by calibration.

## Installation

Python 3.10 or later is recommended.

```bash
git clone https://github.com/zxcvbnm22003Hyl/magnetloss-calculation.git
cd magnetloss-calculation
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e .[dev]
pytest -q
```

## Command-line usage

### Magnetic field

```bash
apexloss field --out turn_fields.csv
```

### Single-strand finite penetration

Using the bundled tabulated material model:

```bash
apexloss strand \
  --Bpk 20.619192 \
  --temperature 4.2 \
  --nr 160 \
  --dt-us 1
```

Using the legacy 4.2 K Kohler/RRR law used for the frozen validation case:

```bash
apexloss strand \
  --Bpk 20.619192 \
  --nr 160 \
  --dt-us 1 \
  --legacy-kohler
```

### Explicit 304-strand FEM

```bash
apexloss secondary-fem \
  --Bpk 20.619192 \
  --temperature 4.2 \
  --h-mm 0.03 \
  --dt-ms 0.25
```

### Whole-magnet fixed-temperature loss

```bash
apexloss whole-magnet \
  --mode fixed \
  --temperature 4.2 \
  --target-order 2 \
  --outdir output_fixed
```

### Whole-magnet adiabatic thermal mapping

```bash
apexloss whole-magnet \
  --mode thermal \
  --temperature 4.2 \
  --target-order 2 \
  --dt-us 20 \
  --outdir output_thermal
```

## Python API

Magnetic field:

```python
from apexloss import MagnetGeometry, FiniteTurnBiotSavart

field = FiniteTurnBiotSavart(
    MagnetGeometry(),
    current_A=47.29e3,
    source_order=20,
)

print(field.center_field_T())
turns = field.turn_center_fields()
```

Single-strand magnetic diffusion:

```python
from apexloss import SingleStrandDiffusion
from apexloss.materials import LegacyKohlerRRR

model = SingleStrandDiffusion(
    diameter_m=0.2e-3,
    radial_cells=160,
)

result = model.simulate_trapezoid(
    B_peak_T=20.619192,
    dt_s=1e-6,
    tail_s=23e-3,
    resistivity=LegacyKohlerRRR(),
)

print(result.Q_total_J_per_m)
```

Explicit 304-strand FEM:

```python
from apexloss.multifilament_fem import MultifilamentFEM

fem = MultifilamentFEM.secondary304(
    mesh_h_m=0.03e-3
)

result = fem.run_trapezoid(
    B_peak_T=20.619192,
    temperature_K=4.2,
    dt_s=0.25e-3,
)

print(result.energy_J_per_m)
```

## Repository layout

```text
src/apexloss/
    field.py                 finite rectangular-turn Biot–Savart
    strand_diffusion.py      single-strand finite penetration
    multifilament_fem.py     explicit 1/19/304-strand 2-D FEM
    surrogate.py             304-strand loss surrogate
    whole_magnet.py          162-turn field-to-loss mapping
    materials.py             rho(T,B), heat capacity, enthalpy
    geometry.py              reference magnet geometry
    waveforms.py             pulse definitions
    data/                    bundled material and surrogate data

examples/                    executable examples
legacy/                      original Python models retained for provenance
reference_results/           frozen regression results
tests/                       unit and regression tests
configs/                     reference configuration
docs/                        numerical-method and validation notes
```

## Reproducibility

Run the compact validation script:

```bash
python scripts/validate_reference.py
```

The repository also includes GitHub Actions CI for supported Python versions.

## Main assumptions and limitations

- strands are electrically insulated;
- inter-strand and inter-subcable coupling-current loss is not included;
- real twist-pitch geometry and helicoidal transformation are not included in v0.1.0;
- structural-metal eddy-current loss is not included;
- joints and leads are not included;
- the thermal model is adiabatic for one pulse;
- the high-purity-Al resistivity model is a design model and should eventually be replaced by measured finished-strand/cable data;
- the 304-strand FEM currently uses an embedded, non-body-fitted mesh;
- the fast thermal surrogate is reduced-order and does not retain the complete local electromagnetic state history.

See:

- [Numerical methods](docs/numerical_methods.md)
- [Model assumptions](docs/model_assumptions.md)
- [Validation snapshot](docs/validation.md)
- [中文说明](docs/README_zh-CN.md)

## Citation

If this code is used in academic work, cite the repository release and the associated methodology paper when available.

## License

BSD-3-Clause.

Copyright © 2026 hery rainforeset.
