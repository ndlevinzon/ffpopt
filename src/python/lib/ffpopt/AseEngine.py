"""Persistent ASE calculators and in-process geomeTRIC (same optimizer settings).

Reloading an ML potential in a fresh ``python -m geometric`` process is the
dominant cost of a constrained scan. This module reuses the unrestrained
base already stored on ``ListOfStruct.calc`` by ``BuildCalc`` and drives
geomeTRIC's ``EngineASE`` in the current process. Convergence keys, coordinate
system, and maxiter are the caller's existing ``GetStandardOptions`` values --
this is not a new optimizer.
"""

from __future__ import annotations

import os
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np

PathLike = Union[str, Path]


def use_geometric_subprocess() -> bool:
    """True when ``FFPOPT_GEOMETRIC_SUBPROCESS=1`` forces the legacy CLI path."""
    raw = os.environ.get("FFPOPT_GEOMETRIC_SUBPROCESS", "")
    return str(raw).strip().lower() not in {"", "0", "false", "no", "off"}


def calc_cache_key(los, struct) -> tuple:
    """Identity for wavefront pool reuse (model / charge / parm).

    Not a second calculator store: the reusable base lives on ``los.calc``.
    Restraint values are omitted so a pool can be reused across scan angles.
    """
    args = getattr(los, "args", None)
    if isinstance(args, Mapping):
        model = str(args.get("model", "sander")).upper()
    else:
        model = str(getattr(args, "model", "sander")).upper()
    charge = None
    parm = None
    try:
        charge = struct.GetCharge()
    except Exception:
        pass
    try:
        parm = struct.data.get("parm")
    except Exception:
        pass
    return (model, charge, parm)


def _wrap_restrained(base, reslist):
    from . ase.calculator import RestrainedCalculator

    return RestrainedCalculator(base, reslist)


def get_persistent_calc(los, struct, reslist=None):
    """Return a calculator using ``los.calc`` as the unrestrained base cache.

    ``ListOfStruct.BuildCalc`` already stores the reusable, unrestrained
    calculator on ``los.calc`` (parm / spin / charge). Restraint wraps are
    returned for this call only and are not written back to ``los.calc``,
    matching ``BuildRestrainedCalc``.
    """
    base = los.BuildCalc(struct)
    if reslist is not None:
        return _wrap_restrained(base, reslist)
    return base


def _normalize_converge(converge) -> Optional[list]:
    if converge is None:
        return None
    if isinstance(converge, str):
        parts = converge.split()
        return parts or None
    if isinstance(converge, Sequence) and not isinstance(converge, (bytes, bytearray)):
        return list(converge)
    return [str(converge)]


def write_plain_xyz(path: PathLike, atoms) -> None:
    """Write element + xyz only (no charge columns) for geomeTRIC Molecule."""
    import ase.io

    ase.io.write(str(path), atoms, format="xyz", parallel=False)


def _rm_path(path: Path) -> bool:
    try:
        if path.is_symlink() or path.is_file():
            path.unlink()
            return True
        if path.is_dir():
            shutil.rmtree(path)
            return True
    except OSError:
        return False
    return False


def prepare_geometric_tmpdir(prefix: PathLike) -> str:
    """Remove leftover ``{prefix}.tmp`` so geomeTRIC ``os.makedirs`` can succeed.

    geomeTRIC calls ``os.makedirs(prefix + '.tmp')`` without ``exist_ok``. A
    killed worker or job restart leaves that directory and raises
    ``FileExistsError``. Also drops ``{prefix}.r*.tmp`` recovery dirs for the
    same prefix.
    """
    prefix = str(prefix)
    tmpdir = prefix + ".tmp"
    _rm_path(Path(tmpdir))
    parent = Path(prefix).parent
    stem = Path(prefix).name
    if parent.is_dir() and stem:
        try:
            for p in parent.iterdir():
                name = p.name
                if name.startswith(stem + ".r") and name.endswith(".tmp"):
                    _rm_path(p)
        except OSError:
            pass
    return tmpdir


def cleanup_geometric_scratch(prefix: PathLike, *, keep_optim: bool = False) -> int:
    """Remove geomeTRIC sidecars for one opt prefix."""
    prefix = os.path.normpath(str(prefix))
    n = 0
    if _rm_path(Path(prefix + ".tmp")):
        n += 1
    for suf in (".log", ".nsf", ".xyz", ".json", ".cons.inp"):
        if _rm_path(Path(prefix + suf)):
            n += 1
    if not keep_optim and _rm_path(Path(prefix + "_optim.xyz")):
        n += 1
    parent = Path(prefix).parent
    stem = Path(prefix).name
    if parent.is_dir():
        try:
            for p in parent.iterdir():
                name = p.name
                if name.startswith(stem + ".r"):
                    _rm_path(p)
                    n += 1
        except OSError:
            pass
    return n


def patch_geometric_tmp_makedirs() -> None:
    """geomeTRIC ``os.makedirs(prefix+'.tmp')`` has no ``exist_ok`` (FileExistsError)."""
    try:
        import geometric.optimize as go
    except ImportError:
        return
    os_mod = getattr(go, "os", None)
    if os_mod is None:
        return
    current = os_mod.makedirs
    if getattr(current, "_ffpopt_exist_ok", False):
        return

    def _makedirs(name, mode=0o777, exist_ok=False):
        text = str(name).rstrip("/\\")
        if text.endswith(".tmp"):
            exist_ok = True
        return current(name, mode=mode, exist_ok=exist_ok)

    _makedirs._ffpopt_exist_ok = True  # type: ignore[attr-defined]
    os_mod.makedirs = _makedirs


def run_geometric_inprocess(
    atoms,
    calc,
    *,
    prefix: PathLike,
    constraints_path: Optional[PathLike] = None,
    coordsys: str = "tric",
    maxiter: int = 500,
    converge="set GAU",
    enforce: Optional[float] = None,
    log_ini: Optional[PathLike] = None,
    **extra_kwargs: Any,
):
    """Run geomeTRIC in-process with an existing ASE calculator.

    Parameters
    ----------
    atoms : ase.Atoms
        Starting geometry (Angstrom).
    calc
        ASE calculator already constructed for this molecule.
    prefix : path-like
        Basename for geometric's log / tmpdir / ``_optim.xyz`` outputs.
    constraints_path : path-like, optional
        Geometric constraint file (same format as the CLI path).
    coordsys, maxiter, converge, enforce
        Forwarded to :func:`geometric.optimize.run_optimizer`.
    log_ini : path-like, optional
        Geometric logging INI.

    Returns
    -------
    dict
        ``coords`` (Ang, ndarray), ``energy_ha`` (Hartree or None),
        ``progress`` (geometric Molecule trajectory).
    """
    patch_geometric_tmp_makedirs()

    from geometric.ase_engine import EngineASE
    from geometric.molecule import Molecule
    from geometric.optimize import run_optimizer

    prefix = str(prefix)
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    xyz_path = prefix + ".xyz"
    write_plain_xyz(xyz_path, atoms)

    M = Molecule(xyz_path)[0]
    engine = EngineASE(M, calc)

    kwargs: dict[str, Any] = {
        "customengine": engine,
        "input": xyz_path,
        "prefix": prefix,
        "coordsys": coordsys,
        "maxiter": int(maxiter),
    }
    conv = _normalize_converge(converge)
    if conv is not None:
        kwargs["converge"] = conv
    if enforce is not None:
        kwargs["enforce"] = float(enforce)
    if constraints_path is not None:
        kwargs["constraints"] = str(constraints_path)
    if log_ini is not None and str(log_ini):
        kwargs["logIni"] = str(log_ini)
    kwargs.update(extra_kwargs)

    prepare_geometric_tmpdir(prefix)
    try:
        progress = run_optimizer(**kwargs)
    except FileExistsError:
        prepare_geometric_tmpdir(prefix)
        progress = run_optimizer(**kwargs)

    coords = np.asarray(progress.xyzs[-1], dtype=float)
    energy_ha = None
    qm_e = getattr(progress, "qm_energies", None)
    if qm_e is not None and len(qm_e) > 0:
        energy_ha = float(qm_e[-1])

    return {
        "coords": coords,
        "energy_ha": energy_ha,
        "progress": progress,
    }
