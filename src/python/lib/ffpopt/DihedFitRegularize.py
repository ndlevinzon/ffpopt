"""Regularized Fourier torsion fits (Tikhonov / truncated SVD + physical caps).

The isolated-torsion linear problem is

    y(φ) = E_HL'(φ) − E_MM_without_this_quartet'(φ)
    V(φ) = Σ_n PK_n (1 + cos(n φ − γ_n))     # Amber PK = V_n / 2

Naive ``pinv`` on periods 1..3 with clustered or noisy leftover invents
large cancelling harmonics. Those then blow up the next relaxed MM scan
(11 → 25 → 82 kcal/mol on CHAPS). This module is the numerically stable
replacement: truncated SVD / Tikhonov, an absolute PK clip, optional AIC
model-order selection, a peak-to-peak barrier cap, and chemical rotor
caps (alkane / sulfate / ammonium) so physically stiff leftover is not
dumped into the scanned torsion.
"""

from __future__ import annotations

import os
from itertools import product
from typing import Any

import numpy as np


def _envf(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw is None or str(raw).strip() == "":
        return float(default)
    return float(raw)


def _env_flag(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


def dihed_fc_abs_max() -> float:
    """Absolute Amber PK bound (kcal/mol). Default 25."""

    return _envf("FFPOPT_DIHED_FC_MAX", 25.0)


def clip_dihed_fcs(fcs, *, where: str = "") -> np.ndarray:
    """Clip Fourier PKs to ±:func:`dihed_fc_abs_max`."""

    cap = dihed_fc_abs_max()
    arr = np.asarray(fcs, dtype=float).copy()
    out = np.clip(arr, -cap, cap)
    if where and np.any(out != arr):
        print(
            f"[fit] clip PKs at {where}: "
            f"{arr.tolist()} -> {out.tolist()} (cap={cap})"
        )
    return out


def tikhonov_svd_solve(
    A: np.ndarray,
    y: np.ndarray,
    lam: float = 0.0,
    rel_cutoff: float = 1.0e-4,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Solve ``min ||A x − y||² + λ||x||²`` with truncated SVD.

    Modes with ``s_i < rel_cutoff * s_max`` are dropped. ``lam=0`` is
    truncated Moore–Penrose (the usual fix for cancelling harmonics).
    """

    A = np.asarray(A, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    if A.ndim != 2:
        raise ValueError(f"A must be 2-D, got {A.shape}")
    if A.shape[0] != y.shape[0]:
        raise ValueError(f"A {A.shape} vs y {y.shape}")
    if A.size == 0:
        return np.zeros(0), {"n_kept": 0, "s": np.zeros(0)}

    u, s, vt = np.linalg.svd(A, full_matrices=False)
    smax = float(s[0]) if s.size else 1.0
    keep = s >= (float(rel_cutoff) * smax) if smax > 0 else np.zeros(s.shape, dtype=bool)
    filt = np.zeros_like(s)
    nz = keep & (s > 0)
    filt[nz] = s[nz] / (s[nz] ** 2 + float(lam))
    x = (vt.T * filt) @ (u.T @ y)
    return x, {
        "n_kept": int(np.sum(keep)),
        "s": s,
        "lam": float(lam),
        "rel_cutoff": float(rel_cutoff),
    }


def dense_torsion_ptp(dfcn, n: int = 361) -> float:
    """Peak-to-peak of ``V(φ)`` on a 1° grid (kcal/mol)."""

    angs = np.linspace(0.0, 360.0, int(n), endpoint=True)
    v = np.asarray(dfcn.CptEne(angs), dtype=float)
    return float(np.max(v) - np.min(v))


def _scale_fcs_to_ptp(dfcn, target: float) -> None:
    ptp = dense_torsion_ptp(dfcn)
    if ptp <= 1.0e-12 or target <= 0:
        return
    scale = float(target) / ptp
    dfcn.SetFCs([float(p.fc) * scale for p in dfcn.prims])


def solve_regularized_fcs(
    A: np.ndarray,
    y: np.ndarray,
    *,
    dfcn=None,
    where: str = "",
    n_fc: int | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Tikhonov solve, clip PKs, optionally scale ``dfcn`` to a barrier cap.

    ``A`` may include a trailing constant column; only the first ``n_fc``
    (or ``len(dfcn.prims)``) entries are treated as Fourier PKs.
    """

    lam = _envf("FFPOPT_DIHED_RIDGE_LAMBDA", 0.0)
    x, info = tikhonov_svd_solve(A, y, lam=lam, rel_cutoff=1.0e-4)
    if n_fc is None:
        n_fc = 0 if dfcn is None else len(dfcn.prims)
    n_fc = max(0, min(int(n_fc), x.size))
    x = x.copy()
    if n_fc:
        x[:n_fc] = clip_dihed_fcs(x[:n_fc], where=where or "regularized")
    if dfcn is not None and n_fc:
        dfcn.SetFCs(x[:n_fc])
        barrier_abs = _envf("FFPOPT_DIHED_BARRIER_ABS", 30.0)
        alpha = _envf("FFPOPT_DIHED_BARRIER_ALPHA", 1.0)
        ptp = dense_torsion_ptp(dfcn)
        info["dense_ptp"] = ptp
        cap = barrier_abs * alpha if alpha > 0 else barrier_abs
        if cap > 0 and ptp > cap:
            _scale_fcs_to_ptp(dfcn, cap)
            x[:n_fc] = np.array([p.fc for p in dfcn.prims], dtype=float)
            info["dense_ptp"] = dense_torsion_ptp(dfcn)
            if where:
                print(
                    f"[fit] barrier cap at {where}: ptp {ptp:.2f} -> "
                    f"{info['dense_ptp']:.2f} kcal/mol (cap={cap})"
                )
        x[:n_fc] = clip_dihed_fcs(x[:n_fc], where=where or "post-barrier")
        dfcn.SetFCs(x[:n_fc])
        info["dense_ptp"] = dense_torsion_ptp(dfcn)
    return x, info


def parse_dihed_type_key(key: str) -> tuple[str, str, str, str]:
    """Split ``{res}_{t1-t2-t3-t4}`` or ``t1-t2-t3-t4`` into four Amber types."""

    s = str(key).strip()
    if "_" in s:
        _left, right = s.split("_", 1)
        if "-" in right:
            s = right
    parts = s.split("-")
    if len(parts) != 4:
        raise ValueError(f"expected 4 Amber types in dihedral key {key!r}")
    return tuple(p.strip() for p in parts)  # type: ignore[return-value]


_ALKANE = frozenset({"c3", "c6", "hc", "h1", "h2", "h3", "ha"})
_SULF_P = frozenset({"s4", "s6", "p", "p4", "p5"})
_AMMON = frozenset({"n4", "n3"})
_ETHER = frozenset({"oh", "os"})
_POLAR_S = frozenset({"ss", "sh"})
_UNSAT = frozenset({"c", "c2", "ca", "n", "ns", "o", "cd", "cc", "ne", "nb"})


def classify_dihed_rotor(type_key: str) -> str:
    """Chemical class of a proper torsion from its Amber type string."""

    t1, t2, t3, t4 = parse_dihed_type_key(type_key)
    types = (t1.lower(), t2.lower(), t3.lower(), t4.lower())
    central = (types[1], types[2])
    if any(t in _SULF_P for t in types):
        return "sulfate_phosphate"
    if any(t in _AMMON for t in types):
        return "amine_ammonium"
    if any(t in _ETHER for t in central):
        return "alcohol_ether"
    if any(t in _POLAR_S for t in types):
        return "polar_sp3"
    if all(t in _ALKANE for t in types):
        return "alkane"
    if any(t in _UNSAT for t in types):
        return "unsaturated"
    return "sp3_sp3"


def apply_sp3_rotor_policy(dfcn, type_key: str, *, where: str = ""):
    """Zero or cap PKs that would give an unphysical isolated-torsion barrier."""

    ptp = dense_torsion_ptp(dfcn)
    kind = classify_dihed_rotor(type_key)
    sp3_max = _envf("FFPOPT_DIHED_SP3_BARRIER_MAX", 20.0)
    alk_max = _envf("FFPOPT_DIHED_ALKANE_BARRIER_MAX", 5.0)
    polar_max = _envf("FFPOPT_DIHED_POLAR_SP3_BARRIER_MAX", 8.0)
    sul_max = _envf("FFPOPT_DIHED_SULFATE_BARRIER_MAX", 10.0)
    sul_cap = _envf("FFPOPT_DIHED_SULFATE_BARRIER_CAP", 4.0)

    action = "keep"
    if kind == "sulfate_phosphate":
        if ptp > sul_max:
            dfcn.SetFCs([0.0] * len(dfcn.prims))
            action = "zero_sulfate_phosphate"
        elif ptp > sul_cap:
            _scale_fcs_to_ptp(dfcn, sul_cap)
            action = "cap_sulfate_phosphate"
    elif kind == "alkane":
        if ptp > sp3_max:
            dfcn.SetFCs([0.0] * len(dfcn.prims))
            action = "zero_alkane"
        elif ptp > alk_max:
            _scale_fcs_to_ptp(dfcn, alk_max)
            action = "cap_alkane"
    elif kind == "amine_ammonium":
        if ptp > polar_max:
            _scale_fcs_to_ptp(dfcn, polar_max)
            action = "cap_amine_ammonium"
    elif kind == "alcohol_ether":
        if ptp > polar_max:
            _scale_fcs_to_ptp(dfcn, polar_max)
            action = "cap_alcohol_ether"
    elif kind == "polar_sp3":
        if ptp > polar_max:
            _scale_fcs_to_ptp(dfcn, polar_max)
            action = "cap_polar_sp3"
    elif kind == "sp3_sp3":
        if ptp > sp3_max:
            dfcn.SetFCs([0.0] * len(dfcn.prims))
            action = "zero_sp3_sp3"

    if where and action != "keep":
        print(
            f"[fit] rotor policy {action} at {where} ({type_key}): "
            f"ptp {ptp:.2f} -> {dense_torsion_ptp(dfcn):.2f} kcal/mol"
        )
    return dfcn, action, ptp


apply_chemical_rotor_policy = apply_sp3_rotor_policy


def nprim_select_enabled() -> bool:
    return _env_flag("FFPOPT_DIHED_NPRIM_SELECT", default=True)


def _aic(rss: float, npts: int, k: int) -> float:
    if npts <= 0:
        return 1.0e30
    return float(npts) * np.log(max(rss, 1.0e-16) / float(npts)) + 2.0 * float(k)


def fit_fourier_nprim(angs, y, nprim_max: int, idxs, pname: str = ""):
    """AIC-select periods ``1..n`` (phase 0) up to ``nprim_max``."""

    from .Dihedrals import GetDihedClasses, MultiDihedFcn, PrimDihedFcn

    angs = np.asarray(angs, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    nprim_max = max(1, int(nprim_max))
    best = None
    for nprim in range(1, nprim_max + 1):
        dfcn = GetDihedClasses(idxs=list(idxs))[nprim][0]
        npts = y.size
        A = np.zeros((npts, nprim + 1))
        for i, prim in enumerate(dfcn.prims):
            A[:, i] = prim.CptEterm(angs)
        A[:, nprim] = 1.0
        x, info = solve_regularized_fcs(
            A, y, dfcn=dfcn, where=f"{pname}:n{nprim}", n_fc=nprim
        )
        v = np.asarray(dfcn.CptEne(angs), dtype=float) + float(x[-1])
        rss = float(np.dot(y - v, y - v))
        aic = _aic(rss, npts, nprim + 1)
        rec = (aic, rss, nprim, dfcn, x, info)
        if best is None or aic < best[0]:
            best = rec
    aic, rss, nprim, dfcn, x, info = best
    info = dict(info)
    info["nprim"] = nprim
    info["aic"] = aic
    info["rss"] = rss
    if pname:
        print(
            f"[fit] AIC nprim={nprim}/{nprim_max} at {pname}: "
            f"rss={rss:.4g} aic={aic:.3f} PKs={[round(p.fc, 4) for p in dfcn.prims]}"
        )
    return dfcn, x, info


def phase_variant_functions(idxs, nprim: int):
    """Amber discrete phases 0°/180° for periods 1..nprim (capped at 16)."""

    from .Dihedrals import MultiDihedFcn, PrimDihedFcn

    nprim = max(1, int(nprim))
    variants = []
    # 2^nprim grows fast; keep the search cheap for nprim>4.
    if nprim > 4:
        variants.append(
            MultiDihedFcn(
                idxs, [PrimDihedFcn(1, 0, n) for n in range(1, nprim + 1)]
            )
        )
        variants.append(
            MultiDihedFcn(
                idxs, [PrimDihedFcn(1, 180, n) for n in range(1, nprim + 1)]
            )
        )
        return variants
    for phases in product((0.0, 180.0), repeat=nprim):
        prims = [PrimDihedFcn(1, ph, n) for n, ph in enumerate(phases, start=1)]
        variants.append(MultiDihedFcn(idxs, prims))
    return variants
