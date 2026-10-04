# APEX-PulseLoss

**APEX-PulseLoss** is a pure-Python research code for magnetic-field calculation and AC-loss estimation in fast pulsed solenoids made from electrically insulated, high-purity metallic multifilament conductors.

The repository contains four independent model levels:

1. **Finite-cross-section Biot–Savart magnet model** for the 162-turn APEX Phase-I solenoid.
2. **Single-strand finite magnetic penetration model** for a round metallic strand in transverse pulsed field.
3. **Explicit 2-D P1 multifilament FEM** for 1, 19, or 304 electrically insulated strands.
4. **304-strand loss surrogate + whole-magnet mapping** with per-turn Gauss sampling and optional adiabatic temperature feedback.

No MATLAB, COMSOL, GetDP, or proprietary solver is required.

> Status: research prototype / reproducible numerical model. It is not a certified engineering or safety-analysis package.

## Physical scope

The default reference configuration is the APEX Phase-I pulsed solenoid:

- 9 double pancakes = 18 single pancakes
- 9 radial turns per single pancake
- 162 total turns
- `Ri = 125 mm`, `Ro = 215.7 mm`
- active axial length `332.7 mm`
- cable envelope `8.3 mm x 16.9 mm` (radial x axial)
- peak current `47.29 kA`
- 3 ms rise / 1 ms flat / 3 ms decay trapezoid
- 2432 x 0.2 mm electrically insulated high-purity Al strands per turn
- bundled material table based on an RRR=3000 design model

The explicit 304-strand local model represents one secondary subcable:

```text
0.2 mm strand -> 19-strand primary -> 16 primaries -> 304-strand secondary
                                                   -> x8 -> 2432-strand turn cable
```

Inter-strand electrical coupling loss is **not** included in the default models because strands are assumed electrically insulated.

## Model hierarchy

### 1. Finite rectangular-turn Biot–Savart model

Each macro turn is treated as a rectangular current-density domain and integrated as a continuum of circular loops with tensor-product Gauss–Legendre quadrature.

This avoids the self-field singularity that occurs when the field is evaluated at the center of a zero-thickness filamentary turn.

Reference result at `47.29 kA`:

```text
B0 ≈ 20.1952 T
maximum turn-center |B| ≈ 20.6301 T
```

### 2. Single-strand finite magnetic penetration

`apexloss.strand_diffusion.SingleStrandDiffusion` solves the time-domain magnetic diffusion problem for the transverse `m=1` mode of a round strand.

It is the refactored version of the original APEX single-strand penetration code and is intentionally kept independent of the 304-strand FEM.

It can be compared directly against the complete-penetration expression

```math
Q_{\rm CP}
=\int \frac{\pi a^4}{4\rho(B,T)}
\left(\frac{dB}{dt}\right)^2dt.
```

The numerical diffusion model retains finite penetration and the post-pulse diffusion tail.

### 3. Explicit 2-D multifilament FEM

`apexloss.multifilament_fem.MultifilamentFEM` uses a scalar `A_z` magnetic-diffusion formulation on a triangular P1 mesh.

For every electrically insulated strand, the strand-average electric field is projected out so that the net transport current is zero in field-only calculations. This isolates strand-internal eddy loss.

Predefined layouts:

- one 0.2 mm strand
- one 19-strand `1+6+12` primary subcable
- one 304-strand `16 x 19` secondary subcable

The current implementation uses an embedded structured triangular mesh. It is intentionally simple and transparent; a body-fitted mesh would reduce the remaining geometric discretization error.

### 4. Whole-magnet surrogate mapping

The bundled database contains the 304-strand full-pulse FEM response on a grid of magnetic field and temperature. Directional dependence is represented as

```math
Q(\theta)=Q_{xx}\cos^2\theta+Q_{yy}\sin^2\theta
+2Q_{xy}\sin\theta\cos\theta.
```

For whole-magnet calculations, each macro turn is sampled using Gauss points over its cross section. The default is `2 x 2 = 4` points per turn, so the 162-turn magnet uses 648 macro field points.

A fixed-temperature whole-magnet loss calculation is therefore

```text
162 finite turns
  -> 648 local (Br, Bz) points
  -> 304-strand surrogate
  -> x8 secondary subcables
  -> integrate over each turn length
  -> whole-magnet loss
```

The package also provides a fast adiabatic thermal mapper derived from the full-pulse surrogate. This is a reduced-order reconstruction and should not be confused with the fully state-resolved 304-strand transient FEM.

## Installation

Clone the repository and install in editable mode:

```bash
git clone <your-repository-url>
cd APEX-PulseLoss
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e .[dev]
pytest -q
```

## Command-line examples

### Magnet field

```bash
apexloss field --out turn_fields.csv
```

### Single-strand finite penetration

Tabulated `rho(T,B)`:

```bash
apexloss strand --Bpk 20.619192 --temperature 4.2 --nr 160 --dt-us 1
```

Legacy 4.2 K Kohler/RRR model used by the original validation script:

```bash
apexloss strand --Bpk 20.619192 --nr 160 --dt-us 1 --legacy-kohler
```

### Explicit 304-strand FEM

```bash
apexloss secondary-fem --Bpk 20.619192 --temperature 4.2 --h-mm 0.03 --dt-ms 0.25
```

This calculation is substantially more expensive than surrogate lookup.

### Whole-magnet fixed-temperature loss

```bash
apexloss whole-magnet --mode fixed --temperature 4.2 --target-order 2 --outdir output_fixed
```

### Whole-magnet adiabatic reduced-order calculation

```bash
apexloss whole-magnet --mode thermal --temperature 4.2 --target-order 2 --dt-us 20 --outdir output_thermal
```

## Python API examples

```python
from apexloss import MagnetGeometry, FiniteTurnBiotSavart

field = FiniteTurnBiotSavart(MagnetGeometry(), current_A=47.29e3, source_order=20)
print(field.center_field_T())
turns = field.turn_center_fields()
```

Single-strand penetration:

```python
from apexloss import SingleStrandDiffusion

model = SingleStrandDiffusion(diameter_m=0.2e-3, radial_cells=160)
result = model.simulate_trapezoid(B_peak_T=20.619192, dt_s=1e-6)
print(result.Q_total_J_per_m)
```

304-strand explicit FEM:

```python
from apexloss.multifilament_fem import MultifilamentFEM

fem = MultifilamentFEM.secondary304(mesh_h_m=0.03e-3)
result = fem.run_trapezoid(B_peak_T=20.619192, temperature_K=4.2, dt_s=0.25e-3)
print(result.energy_J_per_m)
```

## Reference numerical results

The `reference_results/` directory freezes the current validation state. Representative values are:

```text
finite-turn Biot-Savart center field        ~20.1952 T
fixed-4.2 K whole-magnet 304-surrogate loss ~73.4 kJ/pulse (2x2+ Gauss converged)
old 648-point single-strand spatial RVE      18.8466 kJ eddy loss with thermal feedback
304-strand reduced-order thermal model       ~19.9-20.3 kJ eddy loss depending on time-history reconstruction
```

The difference between the independent local models is useful as a model-form uncertainty estimate; it should not be hidden by artificial tuning.

## Important assumptions and limitations

- strands are electrically insulated; no inter-strand or inter-subcable coupling-current loss is included
- no helicoidal transformation or explicit twist-pitch physics is included in v0.1.0
- no structural-metal eddy-current loss
- no joint/lead loss
- adiabatic single-pulse thermal response only
- the bundled compact `rho0(T)` curve plus Kohler reconstruction is a design model; the 20.5 T upper bound is especially important when interpreting higher-field points
- the 304-strand FEM uses a non-body-fitted embedded mesh
- the fast thermal surrogate reconstructs an instantaneous loss coefficient from full-pulse data and is therefore reduced-order, not a replacement for a fully state-resolved local transient FEM

See `docs/model_assumptions.md` for details.

## Repository layout

```text
src/apexloss/
  field.py                 finite rectangular-turn Biot-Savart
  strand_diffusion.py      single-strand finite penetration
  multifilament_fem.py     explicit 1/19/304-strand 2-D FEM
  surrogate.py             304-strand full-pulse surrogate
  whole_magnet.py          162-turn field-to-loss mapping
  materials.py             rho(T,B), heat capacity and enthalpy
  geometry.py              magnet geometry
  waveforms.py             pulse definitions
  data/                    bundled material/surrogate data

examples/                  executable examples
legacy/                    original pure-Python scripts retained verbatim
reference_results/         frozen validation outputs
tests/                     regression/unit tests
configs/                   reference configuration
```

## Reproducibility and provenance

The `legacy/` directory contains the original pure-Python single-strand magnetic-diffusion and full-magnet spatial-RVE scripts from which parts of the refactored package were derived. They are retained for result traceability.

The refactored code intentionally separates:

- macro magnetic field computation
- local strand/multifilament physics
- material constitutive data
- thermal feedback
- whole-magnet quadrature/mapping

This makes each level independently testable.

## License

BSD-3-Clause. Update the copyright holder in `LICENSE` and the author metadata in `CITATION.cff` before a public release if needed.
