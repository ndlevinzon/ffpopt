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
    s_kept = s[keep]
    smin = float(np.min(s_kept)) if s_kept.size else 0.0
    cond = (smax / smin) if smin > 0 else float("inf")
    return x, {
        "n_kept": int(np.sum(keep)),
        "s": s,
        "lam": float(lam),
        "rel_cutoff": float(rel_cutoff),
        "cond": float(cond),
        "smax": smax,
        "smin": smin,
    }


def dense_torsion_ptp(dfcn, n: int = 361) -> float:
    """Peak-to-peak of ``V(φ)`` on a 1° grid (kcal/mol)."""

    angs = np.linspace(0.0, 360.0, int(n), endpoint=True)
    v = np.asarray(dfcn.CptEne(angs), dtype=float)
    return float(np.max(v) - np.min(v))


def scale_fcs_to_ptp(dfcn, target: float) -> None:
    """Uniformly scale PKs so ``V(φ)`` peak-to-peak equals ``target``."""

    ptp = dense_torsion_ptp(dfcn)
    if ptp <= 1.0e-12 or target <= 0:
        return
    scale = float(target) / ptp
    dfcn.SetFCs([float(p.fc) * scale for p in dfcn.prims])


def _scale_fcs_to_ptp(dfcn, target: float) -> None:
    scale_fcs_to_ptp(dfcn, target)


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
    """Cap PKs that would give an unphysical isolated-torsion barrier.

    Never writes all-zero PKs: ``deleteDihedral`` + ``addDihedral(PK=0)``
    removes the original GAFF cosine and leaves a flat MM DIHE panel.
    """

    ptp = dense_torsion_ptp(dfcn)
    kind = classify_dihed_rotor(type_key)
    sp3_max = _envf("FFPOPT_DIHED_SP3_BARRIER_MAX", 20.0)
    alk_max = _envf("FFPOPT_DIHED_ALKANE_BARRIER_MAX", 5.0)
    polar_max = _envf("FFPOPT_DIHED_POLAR_SP3_BARRIER_MAX", 8.0)
    sul_cap = _envf("FFPOPT_DIHED_SULFATE_BARRIER_CAP", 4.0)

    action = "keep"
    if kind == "sulfate_phosphate":
        if ptp > sul_cap:
            scale_fcs_to_ptp(dfcn, sul_cap)
            action = "cap_sulfate_phosphate"
    elif kind == "alkane":
        if ptp > alk_max:
            scale_fcs_to_ptp(dfcn, alk_max)
            action = "cap_alkane"
    elif kind == "amine_ammonium":
        if ptp > polar_max:
            scale_fcs_to_ptp(dfcn, polar_max)
            action = "cap_amine_ammonium"
    elif kind == "alcohol_ether":
        if ptp > polar_max:
            scale_fcs_to_ptp(dfcn, polar_max)
            action = "cap_alcohol_ether"
    elif kind == "polar_sp3":
        if ptp > polar_max:
            scale_fcs_to_ptp(dfcn, polar_max)
            action = "cap_polar_sp3"
    elif kind == "sp3_sp3":
        if ptp > sp3_max:
            scale_fcs_to_ptp(dfcn, sp3_max)
            action = "cap_sp3_sp3"

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
    """AIC-select periods ``1..n`` (signed PK, phase 0) up to ``nprim_max``."""

    return fit_leftover_fourier(angs, y, nprim_max, idxs, pname=pname)


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


def format_prims(dfcn) -> str:
    """Compact ``n=1 PK=... γ=...`` string for logs."""

    if dfcn is None or not getattr(dfcn, "prims", None):
        return "(no Fourier terms)"
    parts = []
    for p in dfcn.prims:
        parts.append(
            f"n={int(p.per)} PK={float(p.fc):.4f} γ={float(p.phase):.0f}"
        )
    return "; ".join(parts)


def design_cosine_matrix(angs, periods, phases=None, instance_angs=None) -> np.ndarray:
    """Design matrix ``[Σ_j cos(n φ_j + γ), 1]``.

    Amber writes ``PK (1 + cos(nφ + γ))``. The extra ``PK`` is a DC term
    collinear with a constant column, so fitting ``1+cos`` plus a constant
    is rank-deficient and invents huge cancelling PKs. Cosine columns
    plus one intercept recover the same shape with a well-posed SVD.

    When ``instance_angs`` is ``(npts, ninst)``, each column is the sum
    over type-equivalent quartets. For a C3 rotor (offsets 0/120/240)
    the n=1 and n=2 columns vanish and only n=3 survives — matching
    the scanned total energy, not one sulfonyl oxygen.
    """

    if instance_angs is None:
        instance_angs = np.asarray(angs, dtype=float).reshape(-1, 1)
    else:
        instance_angs = np.asarray(instance_angs, dtype=float)
        if instance_angs.ndim == 1:
            instance_angs = instance_angs.reshape(-1, 1)
    periods = list(periods)
    n = len(periods)
    A = np.empty((instance_angs.shape[0], n + 1), dtype=float)
    for i, per in enumerate(periods):
        ph = 0.0 if phases is None else float(phases[i])
        A[:, i] = np.sum(
            np.cos(np.deg2rad(float(per) * instance_angs + ph)), axis=1
        )
    A[:, n] = 1.0
    return A


def effective_fourier(dfcn, instance_angs) -> np.ndarray:
    """Sum of ``V(φ_j)`` over type-equivalent instances (kcal/mol)."""

    instance_angs = np.asarray(instance_angs, dtype=float)
    if instance_angs.ndim == 1:
        instance_angs = instance_angs.reshape(-1, 1)
    npts = instance_angs.shape[0]
    if dfcn is None or not getattr(dfcn, "prims", None):
        return np.zeros(npts, dtype=float)
    v = np.zeros(npts, dtype=float)
    for j in range(instance_angs.shape[1]):
        v = v + np.asarray(dfcn.CptEne(instance_angs[:, j]), dtype=float).reshape(-1)
    return v


def cap_effective_ptp(dfcn, instance_angs, target: float, *, where: str = "") -> float:
    """Uniformly scale PKs so instance-sum ``V`` peak-to-peak ≤ ``target``."""

    v = effective_fourier(dfcn, instance_angs)
    ptp = float(np.max(v) - np.min(v)) if v.size else 0.0
    if ptp <= 1.0e-12 or target <= 0 or ptp <= target:
        return ptp
    scale = float(target) / ptp
    dfcn.SetFCs([float(p.fc) * scale for p in dfcn.prims])
    new_ptp = float(
        np.max(effective_fourier(dfcn, instance_angs))
        - np.min(effective_fourier(dfcn, instance_angs))
    )
    if where:
        print(
            f"[fit] effective-ptp cap at {where}: {ptp:.2f} -> {new_ptp:.2f} "
            f"kcal/mol (target={target:.2f})"
        )
    return new_ptp


def column_rel_cutoff() -> float:
    """Drop instance-sum harmonics weaker than this fraction of the strongest."""

    return _envf("FFPOPT_DIHED_COL_REL", 0.05)


def surviving_periods(angs, nprim_max: int, instance_angs=None, rel: float | None = None):
    """Periods whose instance-sum cosine column is not numerically cancelled.

    For a C3 rotor, ``||Σ_j cos(φ_j)||`` and ``||Σ_j cos(2φ_j)||`` are ~0
    while ``||Σ_j cos(3φ_j)||`` is O(√N). Fitting the cancelled columns
    inverts noise into PK ~ 250, which then clips to ±25 and dominates
    the one-quartet plot without moving the scan barrier.
    """

    nprim_max = max(1, int(nprim_max))
    periods = list(range(1, nprim_max + 1))
    A = design_cosine_matrix(angs, periods, instance_angs=instance_angs)
    norms = np.array(
        [float(np.linalg.norm(A[:, i])) for i in range(len(periods))],
        dtype=float,
    )
    if rel is None:
        rel = column_rel_cutoff()
    smax = float(np.max(norms)) if norms.size else 0.0
    kept = [
        int(p)
        for p, n in zip(periods, norms)
        if smax > 0.0 and float(n) >= float(rel) * smax
    ]
    if not kept and periods:
        kept = [int(periods[int(np.argmax(norms))])]
    return kept, {int(p): float(n) for p, n in zip(periods, norms)}


def chemical_barrier_cap(type_key: str) -> float:
    """Peak-to-peak cap (kcal/mol) for the instance-sum torsion on a bond."""

    kind = classify_dihed_rotor(type_key)
    if kind == "sulfate_phosphate":
        return _envf("FFPOPT_DIHED_SULFATE_BARRIER_CAP", 4.0)
    if kind == "alkane":
        return _envf("FFPOPT_DIHED_ALKANE_BARRIER_MAX", 5.0)
    if kind in {"amine_ammonium", "alcohol_ether", "polar_sp3"}:
        return _envf("FFPOPT_DIHED_POLAR_SP3_BARRIER_MAX", 8.0)
    if kind == "sp3_sp3":
        return _envf("FFPOPT_DIHED_SP3_BARRIER_MAX", 20.0)
    return _envf("FFPOPT_DIHED_BARRIER_ABS", 30.0)


def fourier_rss(angs, y, dfcn, instance_angs=None) -> tuple[float, float, np.ndarray, float]:
    """RSS of ``y`` vs ``V(φ)`` after an optimal constant offset."""

    y = np.asarray(y, dtype=float).reshape(-1)
    if instance_angs is not None:
        v = effective_fourier(dfcn, instance_angs)
    elif dfcn is None or not getattr(dfcn, "prims", None):
        v = np.zeros_like(y)
    else:
        v = np.asarray(dfcn.CptEne(angs), dtype=float).reshape(-1)
    c = float(np.mean(y - v))
    r = y - v - c
    rss = float(np.dot(r, r))
    yc = y - float(np.mean(y))
    tss = float(np.dot(yc, yc))
    r2 = (1.0 - rss / tss) if tss > 1.0e-16 else 1.0
    return rss, c, v, r2


def _ridge_lambda_override() -> float | None:
    raw = os.environ.get("FFPOPT_DIHED_RIDGE_LAMBDA")
    if raw is None or str(raw).strip() == "":
        return None
    return float(raw)


def _gcv_lambdas(s: np.ndarray) -> list[float]:
    smax = float(s[0]) if s.size else 1.0
    base = (smax ** 2) if smax > 0 else 1.0
    return [0.0, 1.0e-4 * base, 1.0e-3 * base, 1.0e-2 * base, 1.0e-1 * base]


def gcv_tikhonov_solve(
    A: np.ndarray,
    y: np.ndarray,
    *,
    rel_cutoff: float = 1.0e-4,
    lam: float | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Truncated SVD; GCV-select ``λ`` unless ``lam`` is given."""

    A = np.asarray(A, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    if lam is not None:
        x, info = tikhonov_svd_solve(A, y, lam=lam, rel_cutoff=rel_cutoff)
        info = dict(info)
        info["gcv"] = None
        return x, info

    u, s, _vt = np.linalg.svd(A, full_matrices=False)
    n = float(y.size)
    best = None
    for lam_try in _gcv_lambdas(s):
        x, info = tikhonov_svd_solve(A, y, lam=lam_try, rel_cutoff=rel_cutoff)
        smax = float(s[0]) if s.size else 1.0
        keep = s >= (float(rel_cutoff) * smax) if smax > 0 else np.zeros(s.shape, dtype=bool)
        f = np.zeros_like(s)
        nz = keep & (s > 0)
        f[nz] = (s[nz] ** 2) / (s[nz] ** 2 + float(lam_try))
        r = y - A @ x
        rss = float(np.dot(r, r))
        n_eff = float(np.sum(f))
        denom = max(n - n_eff, 1.0)
        gcv = n * rss / (denom ** 2)
        rec = (gcv, lam_try, x, info)
        if best is None or gcv < best[0]:
            best = rec
    gcv, lam_try, x, info = best
    info = dict(info)
    info["gcv"] = float(gcv)
    info["lam"] = float(lam_try)
    return x, info


def huber_irls_solve(
    A: np.ndarray,
    y: np.ndarray,
    *,
    rel_cutoff: float = 1.0e-4,
    lam: float | None = None,
    max_iter: int = 8,
    c: float = 1.345,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Huber IRLS around a GCV/Tikhonov start (down-weights steric outliers)."""

    A = np.asarray(A, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    x, info = gcv_tikhonov_solve(A, y, rel_cutoff=rel_cutoff, lam=lam)
    lam_fixed = float(info.get("lam", 0.0))
    n_irls = 0
    for _ in range(max(0, int(max_iter))):
        r = y - A @ x
        med = float(np.median(r))
        mad = float(np.median(np.abs(r - med)))
        scale = 1.4826 * mad
        if scale < 1.0e-10:
            break
        u = r / (c * scale)
        w = np.ones_like(r)
        big = np.abs(u) > 1.0
        w[big] = 1.0 / np.abs(u[big])
        sw = np.sqrt(np.clip(w, 1.0e-8, None))
        x_new, info = gcv_tikhonov_solve(
            A * sw[:, None],
            sw * y,
            rel_cutoff=rel_cutoff,
            lam=lam_fixed,
        )
        n_irls += 1
        if float(np.max(np.abs(x_new - x))) < 1.0e-8:
            x = x_new
            break
        x = x_new
    info = dict(info)
    info["n_irls"] = n_irls
    info["lam"] = lam_fixed
    return x, info


def _solve_leftover_linear(A: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    lam = _ridge_lambda_override()
    if _env_flag("FFPOPT_DIHED_IRLS", default=True):
        return huber_irls_solve(A, y, lam=lam)
    return gcv_tikhonov_solve(A, y, lam=lam)


def fit_leftover_fourier(
    angs, y, nprim_max: int, idxs, pname: str = "", instance_angs=None
):
    """Fit Amber PKs to isolated leftover with cosine columns, GCV, AIC.

    Signed PK at phase 0° is equivalent (up to a constant) to a 180°
    term, so the old ``2^nprim`` phase enumeration is unnecessary.
    Only instance-sum columns that survive :func:`surviving_periods`
    are fitted, so C3 rotors get n=3 rather than clipped ±25 on n=1,2.
    """

    from .Dihedrals import MultiDihedFcn, PrimDihedFcn

    angs = np.asarray(angs, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    nprim_max = max(1, int(nprim_max))
    kept, norms = surviving_periods(angs, nprim_max, instance_angs=instance_angs)
    if pname:
        bits = ", ".join(f"n={p} ||cos||={norms[p]:.3g}" for p in sorted(norms))
        print(f"[fit] instance-sum column norms at {pname}: {bits}")
        dropped = [p for p in sorted(norms) if p not in kept]
        if dropped:
            print(
                f"[fit] dropping cancelled periods {dropped} at {pname} "
                f"(||cos|| < {column_rel_cutoff():.2g}× strongest {max(norms.values()):.3g})"
            )
    if not kept:
        kept = [1]
    if nprim_select_enabled() and len(kept) > 1:
        candidates = [kept[:i] for i in range(1, len(kept) + 1)]
    else:
        candidates = [kept]
    best = None
    for periods in candidates:
        dfcn = MultiDihedFcn(
            list(idxs), [PrimDihedFcn(1.0, 0.0, n) for n in periods]
        )
        A = design_cosine_matrix(angs, periods, instance_angs=instance_angs)
        x, info = _solve_leftover_linear(A, y)
        n_fc = len(dfcn.prims)
        tag = f"{pname}:n{','.join(str(p) for p in periods)}" if pname else ""
        pks = clip_dihed_fcs(x[:n_fc], where=tag)
        dfcn.SetFCs(pks)
        rss, const, _v, r2 = fourier_rss(angs, y, dfcn, instance_angs=instance_angs)
        aic = _aic(rss, y.size, len(periods) + 1)
        rec = (aic, rss, len(periods), dfcn, x, info, const, r2, tuple(periods))
        if best is None or aic < best[0]:
            best = rec
    aic, rss, nprim, dfcn, x, info, const, r2, periods = best
    info = dict(info)
    info["nprim"] = int(nprim)
    info["periods"] = list(periods)
    info["aic"] = float(aic)
    info["rss"] = float(rss)
    info["r2"] = float(r2)
    info["const"] = float(const)
    info["dropped_periods"] = [p for p in sorted(norms) if p not in kept]
    if pname:
        print(
            f"[fit] leftover LS {pname}: periods={list(periods)} "
            f"rss={rss:.4g} r²={r2:.4f} aic={aic:.3f} "
            f"λ={info.get('lam', 0):.3g} cond={info.get('cond', float('nan')):.3g} "
            f"n_irls={info.get('n_irls', 0)} PKs={[round(p.fc, 4) for p in dfcn.prims]}"
        )
        if float(info.get("cond") or 0.0) > 1.0e6:
            print(
                f"[fit] WARNING {pname}: design-matrix cond={info['cond']:.3g} "
                "(clustered φ or rank deficiency); Tikhonov/GCV truncated the null space"
            )
    return dfcn, x, info


def append_fit_trace(rec: dict, path: str = "fit_trace.jsonl") -> None:
    """Append one JSON object to the fragment-local fit audit log."""

    import json

    clean = {}
    for key, val in rec.items():
        if isinstance(val, np.ndarray):
            clean[key] = [round(float(v), 6) for v in np.ravel(val)[:64]]
        else:
            clean[key] = val
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(clean, default=str) + "\n")
