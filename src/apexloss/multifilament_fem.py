from __future__ import annotations

from dataclasses import dataclass
import math
import time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu

from .constants import MU0
from .materials import ResistivityTable
from .waveforms import TrapezoidPulse


def hex19_centers(strand_diameter_m=0.2e-3):
    d = strand_diameter_m
    pts = []
    for q in range(-2, 3):
        for r in range(-2, 3):
            s = -q - r
            if max(abs(q), abs(r), abs(s)) <= 2:
                pts.append((d * (q + 0.5 * r), d * (math.sqrt(3) / 2 * r)))
    pts = np.asarray(pts, float)
    rr = np.hypot(pts[:, 0], pts[:, 1])
    ang = np.arctan2(pts[:, 1], pts[:, 0])
    return pts[np.lexsort((ang, rr))]


def secondary304_centers(
    strand_diameter_m=0.2e-3,
    primary_pitch_x_m=1.05e-3,
    primary_pitch_y_m=1.05e-3,
):
    base = hex19_centers(strand_diameter_m)
    centers = []
    groups = []
    coords = [-1.5, -0.5, 0.5, 1.5]
    gid = 0
    for jy in coords:
        for ix in coords:
            shift = np.array([ix * primary_pitch_x_m, jy * primary_pitch_y_m])
            for p in base:
                centers.append(shift + p)
                groups.append(gid)
            gid += 1
    return np.asarray(centers), np.asarray(groups, int)


def _structured_tri_mesh(xmin, xmax, ymin, ymax, h):
    nx = int(np.ceil((xmax - xmin) / h))
    ny = int(np.ceil((ymax - ymin) / h))
    xs = np.linspace(xmin, xmax, nx + 1)
    ys = np.linspace(ymin, ymax, ny + 1)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    nodes = np.column_stack([X.ravel(), Y.ravel()])
    tris = []
    for j in range(ny):
        i = np.arange(nx)
        n00 = j * (nx + 1) + i
        n10 = n00 + 1
        n01 = n00 + (nx + 1)
        n11 = n01 + 1
        tris.append(np.column_stack([n00, n10, n11]))
        tris.append(np.column_stack([n00, n11, n01]))
    tris = np.vstack(tris).astype(np.int32)
    boundary = np.unique(
        np.concatenate(
            [
                np.arange(nx + 1),
                ny * (nx + 1) + np.arange(nx + 1),
                np.arange(ny + 1) * (nx + 1),
                np.arange(ny + 1) * (nx + 1) + nx,
            ]
        )
    ).astype(np.int32)
    return nodes, tris, boundary


def _assemble(nodes, tris, centers, radius):
    nnode = len(nodes)
    nstr = len(centers)
    p0 = nodes[tris[:, 0]]
    p1 = nodes[tris[:, 1]]
    p2 = nodes[tris[:, 2]]
    det = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0]) * (
        p1[:, 1] - p0[:, 1]
    )
    area = np.abs(det) / 2
    b = np.column_stack(
        [p1[:, 1] - p2[:, 1], p2[:, 1] - p0[:, 1], p0[:, 1] - p1[:, 1]]
    ) / det[:, None]
    c = np.column_stack(
        [p2[:, 0] - p1[:, 0], p0[:, 0] - p2[:, 0], p1[:, 0] - p0[:, 0]]
    ) / det[:, None]
    kval = (area[:, None, None] * (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :])).reshape(-1)
    rows = np.repeat(tris, 3, axis=1).reshape(-1)
    cols = np.tile(tris, (1, 3)).reshape(-1)
    K = coo_matrix((kval, (rows, cols)), shape=(nnode, nnode)).tocsr()

    cent = (p0 + p1 + p2) / 3
    tree = cKDTree(centers)
    dist, sid = tree.query(cent, k=1)
    cond = dist <= radius
    sidc = sid[cond].astype(int)
    tc = tris[cond]
    ac = area[cond]

    Mloc = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], float) / 12.0
    mval = (ac[:, None, None] * Mloc[None, :, :]).reshape(-1)
    mrows = np.repeat(tc, 3, axis=1).reshape(-1)
    mcols = np.tile(tc, (1, 3)).reshape(-1)
    M = coo_matrix((mval, (mrows, mcols)), shape=(nnode, nnode)).tocsr()

    grows = tc.reshape(-1)
    gcols = np.repeat(sidc, 3)
    gdata = np.repeat(ac / 3.0, 3)
    G = coo_matrix((gdata, (grows, gcols)), shape=(nnode, nstr)).tocsr()
    S = np.bincount(sidc, weights=ac, minlength=nstr)
    if np.any(S <= 0):
        raise RuntimeError("Some strands contain no conductor elements; refine the mesh.")
    C = (G @ diags(1.0 / S) @ G.T).tocsr()
    Mp = (M - C).tocsr()
    Mp.eliminate_zeros()
    return K, Mp, (tc, ac, sidc), S


@dataclass
class MultifilamentFEMResult:
    time_s: np.ndarray
    power_W_per_m: np.ndarray
    energy_J_per_m: float
    group_energy_J_per_m: np.ndarray
    n_nodes: int
    n_triangles: int
    runtime_s: float
    mean_strand_area_error: float


class MultifilamentFEM:
    """2-D P1 magnetic-diffusion FEM for electrically insulated strands."""

    def __init__(
        self,
        centers_m,
        groups=None,
        strand_diameter_m=0.2e-3,
        domain_half_m=2.8e-3,
        mesh_h_m=0.03e-3,
    ):
        self.centers = np.asarray(centers_m, float)
        self.groups = np.zeros(len(self.centers), dtype=int) if groups is None else np.asarray(groups, int)
        self.radius = strand_diameter_m / 2
        self.domain_half = domain_half_m
        self.mesh_h = mesh_h_m
        self.nodes, self.tris, self.boundary = _structured_tri_mesh(
            -domain_half_m, domain_half_m, -domain_half_m, domain_half_m, mesh_h_m
        )
        self.K, self.Mp, self._element_data, self._strand_area = _assemble(
            self.nodes, self.tris, self.centers, self.radius
        )
        bmask = np.zeros(len(self.nodes), bool)
        bmask[self.boundary] = True
        self.free = np.flatnonzero(~bmask)
        self.Kff = self.K[self.free][:, self.free].tocsc()
        self.Kfb = self.K[self.free][:, self.boundary].tocsr()
        self.Mff = self.Mp[self.free][:, self.free].tocsc()
        self.Mfb = self.Mp[self.free][:, self.boundary].tocsr()

    @classmethod
    def secondary304(cls, **kwargs):
        strand_diameter_m = kwargs.pop("strand_diameter_m", 0.2e-3)
        centers, groups = secondary304_centers(
            strand_diameter_m,
            kwargs.pop("primary_pitch_x_m", 1.05e-3),
            kwargs.pop("primary_pitch_y_m", 1.05e-3),
        )
        return cls(centers, groups, strand_diameter_m=strand_diameter_m, **kwargs)

    @classmethod
    def primary19(cls, strand_diameter_m=0.2e-3, **kwargs):
        c = hex19_centers(strand_diameter_m)
        return cls(c, np.zeros(19, int), strand_diameter_m=strand_diameter_m, **kwargs)

    @classmethod
    def single(cls, strand_diameter_m=0.2e-3, **kwargs):
        return cls(np.array([[0.0, 0.0]]), np.array([0]), strand_diameter_m=strand_diameter_m, **kwargs)

    def run_trapezoid(
        self,
        B_peak_T=20.619192,
        angle_rad=np.pi / 2,
        temperature_K=4.2,
        pulse=TrapezoidPulse(),
        dt_s=0.25e-3,
        resistivity=None,
    ) -> MultifilamentFEMResult:
        table = ResistivityTable() if resistivity is None else resistivity
        n = len(self.nodes)
        A = np.zeros(n)
        nstep = int(round(pulse.duration_s / dt_s))
        times = np.linspace(0, pulse.duration_s, nstep + 1)
        amp = pulse.amplitude(times)
        Bx = B_peak_T * np.cos(angle_rad) * amp
        By = B_peak_T * np.sin(angle_rad) * amp
        xb = self.nodes[self.boundary, 0]
        yb = self.nodes[self.boundary, 1]
        P = np.zeros(nstep)
        Qg = np.zeros(self.groups.max() + 1)
        tc, ac, sidc = self._element_data
        tri_group = self.groups[sidc]
        cache = {}
        t0 = time.time()

        for k in range(nstep):
            Bm = B_peak_T * pulse.amplitude(0.5 * (times[k] + times[k + 1]))
            rho = float(table(temperature_K, abs(Bm))) if isinstance(table, ResistivityTable) else float(table(abs(Bm)))
            sigma = 1.0 / rho
            ccoef = MU0 * sigma / dt_s
            key = round(ccoef, 4)
            if key not in cache:
                cache[key] = splu((self.Kff + ccoef * self.Mff).tocsc())
            Ab = Bx[k + 1] * yb - By[k + 1] * xb
            rhs_full = ccoef * (self.Mp @ A)
            rhs = rhs_full[self.free] - (self.Kfb + ccoef * self.Mfb) @ Ab
            Anew = np.empty_like(A)
            Anew[self.boundary] = Ab
            Anew[self.free] = cache[key].solve(rhs)
            dA = (Anew - A) / dt_s
            P[k] = sigma * float(dA @ (self.Mp @ dA))

            davg_nodes = dA[tc]
            strand_int = np.bincount(
                sidc, weights=ac / 3.0 * davg_nodes.sum(axis=1), minlength=len(self.centers)
            )
            strand_mean = strand_int / self._strand_area
            fvals = davg_nodes - strand_mean[sidc, None]
            sq = ac / 6.0 * (
                np.sum(fvals * fvals, axis=1)
                + fvals[:, 0] * fvals[:, 1]
                + fvals[:, 1] * fvals[:, 2]
                + fvals[:, 2] * fvals[:, 0]
            )
            Pe = sigma * sq
            Qg += np.bincount(tri_group, weights=Pe, minlength=len(Qg)) * dt_s
            A = Anew

        exact_area = np.pi * self.radius**2
        area_err = float(np.mean(np.abs(self._strand_area - exact_area) / exact_area))
        return MultifilamentFEMResult(
            time_s=0.5 * (times[:-1] + times[1:]),
            power_W_per_m=P,
            energy_J_per_m=float(P.sum() * dt_s),
            group_energy_J_per_m=Qg,
            n_nodes=len(self.nodes),
            n_triangles=len(self.tris),
            runtime_s=time.time() - t0,
            mean_strand_area_error=area_err,
        )
