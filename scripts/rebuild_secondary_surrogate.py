"""Rebuild the bundled 304-strand full-pulse loss surrogate from explicit FEM.

This is intentionally a long-running reproducibility script. For each (T,B) point it
runs 0, 45 and 90 degree field orientations and reconstructs Qxx,Qyy,Qxy.
"""
from __future__ import annotations

import argparse
import time
import numpy as np
import pandas as pd

from apexloss.multifilament_fem import MultifilamentFEM


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="secondary304_surrogate_rebuilt.csv")
    ap.add_argument("--h-mm", type=float, default=0.03)
    ap.add_argument("--dt-ms", type=float, default=0.25)
    args = ap.parse_args()

    Bgrid = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20.5, 22.0]
    Tgrid = [4.2, 20.0, 40.0, 77.0]
    fem = MultifilamentFEM.secondary304(mesh_h_m=args.h_mm * 1e-3)
    rows = []

    for T in Tgrid:
        for B in Bgrid:
            if B == 0:
                q0 = q45 = q90 = 0.0
                runtime = 0.0
            else:
                t0 = time.time()
                q0 = fem.run_trapezoid(B, 0.0, T, dt_s=args.dt_ms * 1e-3).energy_J_per_m
                q45 = fem.run_trapezoid(B, np.pi / 4, T, dt_s=args.dt_ms * 1e-3).energy_J_per_m
                q90 = fem.run_trapezoid(B, np.pi / 2, T, dt_s=args.dt_ms * 1e-3).energy_J_per_m
                runtime = time.time() - t0
            qxy = q45 - 0.5 * (q0 + q90)
            qmean = 0.5 * (q0 + q90)
            spread = 0.0 if qmean == 0 else (max(q0, q45, q90) - min(q0, q45, q90)) / qmean
            rows.append(
                {
                    "T_K": T,
                    "Bpk_T": B,
                    "Qxx_Jpm": q0,
                    "Qyy_Jpm": q90,
                    "Qxy_Jpm": qxy,
                    "Qtheta45_Jpm": q45,
                    "Qmean_Jpm": qmean,
                    "angle_spread_rel": spread,
                    "runtime_s": runtime,
                    "B_extrapolated": B > 20.5,
                }
            )
            print(rows[-1])
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(args.out)


if __name__ == "__main__":
    main()
