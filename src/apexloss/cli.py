from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from .field import FiniteTurnBiotSavart
from .geometry import MagnetGeometry
from .materials import LegacyKohlerRRR, ResistivityTable
from .multifilament_fem import MultifilamentFEM
from .strand_diffusion import SingleStrandDiffusion
from .surrogate import SecondaryLossSurrogate
from .whole_magnet import fixed_temperature_loss, QuasiStaticThermalMapper


def _write_summary(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2), encoding="utf-8")


def cmd_field(args):
    model = FiniteTurnBiotSavart(MagnetGeometry(), current_A=args.current_kA * 1e3, source_order=args.source_order)
    df = model.turn_center_fields()
    df.to_csv(args.out, index=False)
    print(f"B0 = {model.center_field_T():.9f} T")
    print(f"turn-center Bmax = {df.Bmag_pk_T.max():.9f} T")
    print(args.out)


def cmd_strand(args):
    model = SingleStrandDiffusion(diameter_m=args.diameter_mm * 1e-3, radial_cells=args.nr)
    if args.legacy_kohler:
        rho = LegacyKohlerRRR()
    else:
        rho = ResistivityTable()
    result = model.simulate_trapezoid(
        B_peak_T=args.Bpk,
        dt_s=args.dt_us * 1e-6,
        tail_s=args.tail_ms * 1e-3,
        temperature_K=args.temperature,
        resistivity=rho,
    )
    d = result.__dict__
    _write_summary(args.out, d)
    print(json.dumps(d, indent=2))


def cmd_secondary(args):
    fem = MultifilamentFEM.secondary304(mesh_h_m=args.h_mm * 1e-3, domain_half_m=args.domain_half_mm * 1e-3)
    result = fem.run_trapezoid(
        B_peak_T=args.Bpk,
        angle_rad=args.angle_deg * 3.141592653589793 / 180,
        temperature_K=args.temperature,
        dt_s=args.dt_ms * 1e-3,
    )
    d = {
        "energy_J_per_m": result.energy_J_per_m,
        "group_energy_J_per_m": result.group_energy_J_per_m.tolist(),
        "n_nodes": result.n_nodes,
        "n_triangles": result.n_triangles,
        "runtime_s": result.runtime_s,
        "mean_strand_area_error": result.mean_strand_area_error,
    }
    _write_summary(args.out, d)
    print(json.dumps(d, indent=2))


def cmd_whole(args):
    geom = MagnetGeometry()
    field = FiniteTurnBiotSavart(geom, current_A=args.current_kA * 1e3, source_order=args.source_order)
    surrogate = SecondaryLossSurrogate()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if args.mode == "fixed":
        turns, Q = fixed_temperature_loss(field, surrogate, args.temperature, args.target_order)
        turns.to_csv(outdir / "turn_loss.csv", index=False)
        summary = {"Q_eddy_J": Q, "temperature_K": args.temperature, "mode": "fixed"}
    else:
        mapper = QuasiStaticThermalMapper(surrogate=surrogate)
        samples, turns, hist, result = mapper.run(
            field,
            target_order=args.target_order,
            dt_s=args.dt_us * 1e-6,
            initial_temperature_K=args.temperature,
            current_peak_A=args.current_kA * 1e3,
        )
        samples.to_csv(outdir / "local_samples.csv", index=False)
        turns.to_csv(outdir / "turn_loss.csv", index=False)
        hist.to_csv(outdir / "time_history.csv", index=False)
        summary = result.__dict__ | {"mode": "quasistatic-thermal"}
    _write_summary(outdir / "summary.json", summary)
    print(json.dumps(summary, indent=2))


def build_parser():
    p = argparse.ArgumentParser(prog="apexloss", description="Pure-Python APEX pulsed-magnet loss models")
    sub = p.add_subparsers(dest="command", required=True)

    f = sub.add_parser("field", help="Compute 162-turn finite-cross-section Biot-Savart fields")
    f.add_argument("--current-kA", type=float, default=47.29)
    f.add_argument("--source-order", type=int, default=20)
    f.add_argument("--out", default="turn_fields.csv")
    f.set_defaults(func=cmd_field)

    s = sub.add_parser("strand", help="Single-strand finite magnetic penetration/diffusion")
    s.add_argument("--Bpk", type=float, default=20.619192)
    s.add_argument("--temperature", type=float, default=4.2)
    s.add_argument("--diameter-mm", type=float, default=0.2)
    s.add_argument("--nr", type=int, default=160)
    s.add_argument("--dt-us", type=float, default=1.0)
    s.add_argument("--tail-ms", type=float, default=23.0)
    s.add_argument("--legacy-kohler", action="store_true")
    s.add_argument("--out", default="strand_diffusion.json")
    s.set_defaults(func=cmd_strand)

    m = sub.add_parser("secondary-fem", help="Explicit 304-strand 2-D field-only FEM")
    m.add_argument("--Bpk", type=float, default=20.619192)
    m.add_argument("--temperature", type=float, default=4.2)
    m.add_argument("--angle-deg", type=float, default=90.0)
    m.add_argument("--h-mm", type=float, default=0.03)
    m.add_argument("--dt-ms", type=float, default=0.25)
    m.add_argument("--domain-half-mm", type=float, default=2.8)
    m.add_argument("--out", default="secondary304_fem.json")
    m.set_defaults(func=cmd_secondary)

    w = sub.add_parser("whole-magnet", help="Whole-magnet field-to-loss mapping")
    w.add_argument("--mode", choices=["fixed", "thermal"], default="fixed")
    w.add_argument("--temperature", type=float, default=4.2)
    w.add_argument("--current-kA", type=float, default=47.29)
    w.add_argument("--source-order", type=int, default=20)
    w.add_argument("--target-order", type=int, default=2)
    w.add_argument("--dt-us", type=float, default=20.0)
    w.add_argument("--outdir", default="apexloss_output")
    w.set_defaults(func=cmd_whole)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)
