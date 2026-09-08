"""Wavefront IPC, checkpoint, and clash-precheck helpers (same search, less I/O).

Workers receive slim ``{angle, coords}`` jobs and keep ``los`` from the pool
initializer so ``get_persistent_calc`` survives across nodes. Checkpoints strip
live calculators and force arrays. Clash precheck is the original nonbonded
``min_dist`` test, vectorized.
"""

from __future__ import annotations

import os
import pickle
import shutil
from pathlib import Path
from typing import Any, Optional, Sequence, Union

import numpy as np

PathLike = Union[str, Path]


def require_main_guard_for_spawn(api_name: str = "this call") -> None:
    """Abort if a spawn worker re-imported the caller's script."""
    import sys

    depth = 1
    while True:
        try:
            frame = sys._getframe(depth)
        except ValueError:
            return
        if frame.f_globals.get("__name__") == "__mp_main__":
            raise RuntimeError(
                f"{api_name} was re-invoked by a multiprocessing spawn worker "
                "re-importing the calling script. Wrap the call in "
                "`if __name__ == '__main__':` so the worker re-import doesn't "
                "re-execute it. See "
                "https://docs.python.org/3/library/multiprocessing.html"
                "#multiprocessing-programming"
            )
        depth += 1
        if depth > 32:
            return


def wf_checkpoint_every(nproc: int) -> int:
    """How many completed nodes between wavefront checkpoints.

    Default is ``nproc`` (same as before). Override with
    ``FFPOPT_WF_CHECKPOINT_EVERY``.
    """
    raw = os.environ.get("FFPOPT_WF_CHECKPOINT_EVERY")
    if raw is not None and str(raw).strip():
        try:
            return max(1, int(str(raw).strip()))
        except ValueError:
            pass
    return max(1, int(nproc))


def clone_struct_geometry(struct, coords, ene=0.0, frcs=None):
    """Prefer ``Struct.clone_geometry``; fall back to deepcopy for test doubles."""
    import copy

    clone = getattr(struct, "clone_geometry", None)
    if callable(clone):
        return clone(coords=coords, ene=ene, frcs=frcs)
    out = copy.deepcopy(struct)
    updater = getattr(out, "Update", None)
    if callable(updater):
        updater(ene, np.asarray(coords, dtype=float), frcs)
    elif getattr(out, "data", None) is not None:
        out.data["positions"] = np.asarray(coords, dtype=float).tolist()
        out.data["energy"] = ene
    return out


def clear_los_calc(los) -> None:
    """Unbind ``los.calc`` on this copy only (worker init / checkpoint dump).

    Does not ``reset()`` or clear the parent cache. Spawn pickle already
    sends ``calc=None`` via ``ListOfStruct.__getstate__``; the parent keeps
    ``los.calc``. Workers refill it with ``BuildCalc`` on first use.
    """
    if los is None:
        return
    los.calc = None


def _replace_with_retry(src: Path, dst: Path, *, attempts: int = 8) -> None:
    last: Optional[OSError] = None
    for i in range(max(1, int(attempts))):
        try:
            os.replace(src, dst)
            return
        except OSError as exc:
            last = exc
            import time

            time.sleep(0.05 * (2 ** i))
            dst.parent.mkdir(parents=True, exist_ok=True)
    if last is None:
        raise FileNotFoundError(str(src))
    try:
        shutil.copyfile(src, dst)
        src.unlink()
    except OSError:
        raise last from last


def atomic_pickle_dump(obj: Any, path) -> None:
    """Write a pickle via unique ``tmp`` + ``os.replace`` (crash-safe)."""
    import time

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    replaced = False
    try:
        with open(tmp, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
            f.flush()
            os.fsync(f.fileno())
        _replace_with_retry(tmp, path)
        replaced = True
    finally:
        if not replaced:
            try:
                if tmp.exists():
                    tmp.unlink()
            except OSError:
                pass


def write_node_pickle(node: Any, *, verbose: bool = False) -> None:
    """Pickle ``node`` without ``los`` (restored after write)."""
    if verbose:
        print(f"Saving node {node.node_pkl}  (exists? {Path(node.node_pkl).is_file()})")
    los = getattr(node, "los", None)
    node.los = None
    try:
        atomic_pickle_dump(node, node.node_pkl)
    except OSError as exc:
        print(
            f"could not write node pickle {node.node_pkl}: "
            f"{type(exc).__name__}: {exc}"
        )
    finally:
        node.los = los


def pickle_checkpoint_keep_calc_cache(obj: Any, path, los) -> None:
    """Pickle a wavefront after unbinding ``los.calc``, then restore it.

    Checkpoints must not serialize sander/XTB/DFT handles. The parent still
    needs ``los.calc`` so serial ``nproc=1`` does not reload the model after
    every snapshot. Unbind, dump, rebind.
    """
    calc = getattr(los, "calc", None) if los is not None else None
    if los is not None:
        clear_los_calc(los)
    try:
        atomic_pickle_dump(obj, path)
    finally:
        if los is not None and calc is not None:
            los.calc = calc


def slim_node_result(node: Any) -> dict:
    """Build the slim multiprocessing result payload for a completed node."""
    coords = None
    opt = getattr(node, "opt_geom", None)
    data = getattr(opt, "data", None) if opt is not None else None
    if isinstance(data, dict) and data.get("positions") is not None:
        coords = np.asarray(data["positions"], dtype=float)
    return {
        "energy": node.energy,
        "forces": getattr(node, "forces", None),
        "coords": coords,
        "complete": node.complete,
        "error": getattr(node, "error", None),
        "active": node.active,
    }


def apply_slim_node_result(node: Any, result: dict) -> None:
    """Merge a slim worker result into a parent-side node."""
    node.energy = result.get("energy")
    if result.get("forces") is not None:
        node.forces = result["forces"]
    node.complete = bool(result.get("complete", node.complete))
    node.error = result.get("error", getattr(node, "error", None))
    if "active" in result:
        node.active = bool(result["active"])
    coords = result.get("coords")
    if coords is not None:
        node.opt_geom = clone_struct_geometry(
            node.struct, coords, ene=node.energy, frcs=result.get("forces")
        )


def slim_completed_nodes_for_checkpoint(wavefront: Any) -> None:
    """Drop bulky force arrays from completed nodes before pickling."""
    for level in getattr(wavefront, "levels", []) or []:
        for node in getattr(level, "nodes", []) or []:
            if not getattr(node, "complete", False):
                continue
            opt = getattr(node, "opt_geom", None)
            if opt is not None and getattr(opt, "data", None) is not None:
                if "forces" in opt.data:
                    opt.data["forces"] = None
            n_atoms = len(node.struct.data["elements"])
            node.forces = np.zeros((n_atoms, 3))


def replace_node_with_pickle(node: Any, *, found_msg: Optional[str] = None) -> None:
    """Replace node fields from a sidecar pickle if present (restores ``los``)."""
    filename = Path(f"{node.node_pkl}")
    if not filename.is_file():
        return
    los = getattr(node, "los", None)
    msg = found_msg or f"Found existing pickle file for node: {node.node_id}"
    print(msg)
    with open(filename, "rb") as f:
        loaded_node = pickle.load(f)
    node.__dict__.update(loaded_node.__dict__)
    if getattr(node, "los", None) is None:
        node.los = los
    print("Node data replaced with pickle data.")


def has_nonbonded_clash(positions, bonds, min_dist: float = 0.8):
    """Return ``(clashed, i, j, dist)`` for the first nonbonded pair below ``min_dist``.

    Same rule as the original nested ``get_distance`` loop: skip bonded pairs,
    flag any other pair closer than ``min_dist``. Bonded skip is order-independent.
    """
    pos = np.asarray(positions, dtype=float)
    n = int(pos.shape[0])
    if n < 2:
        return False, None, None, None
    bonded = set()
    for b in bonds or []:
        try:
            i, j = int(b[0]), int(b[1])
        except (TypeError, ValueError, IndexError):
            continue
        bonded.add((i, j) if i < j else (j, i))
    diff = pos[:, None, :] - pos[None, :, :]
    dist = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
    iu = np.triu_indices(n, k=1)
    for i, j in zip(iu[0].tolist(), iu[1].tolist()):
        pair = (i, j)
        if pair in bonded:
            continue
        d = float(dist[i, j])
        if d < min_dist:
            return True, i, j, d
    return False, None, None, None


def geometric_prefix_from_node_pkl(node_pkl: PathLike) -> str:
    return str(Path(node_pkl).with_suffix("")) + "_geom"


def sweep_geometric_scratch_dir(
    directory: PathLike,
    *,
    keep_optim_prefixes: Optional[Sequence[PathLike]] = None,
    recursive: bool = False,
) -> int:
    """Remove leftover geomeTRIC ``.nsf`` / ``*.tmp`` / ``*_geom*`` sidecars."""
    from . AseEngine import _rm_path

    directory = Path(directory)
    if not directory.is_dir():
        return 0
    keep = {os.path.normpath(str(p)) for p in (keep_optim_prefixes or [])}
    n = 0

    def _keep_optim_file(path: Path) -> bool:
        if not path.name.endswith("_optim.xyz"):
            return False
        pref = os.path.normpath(str(path)[: -len("_optim.xyz")])
        return pref in keep

    for root, dirs, files in os.walk(directory, topdown=True):
        root_p = Path(root)
        next_dirs = []
        for d in dirs:
            dp = root_p / d
            if d.endswith(".tmp"):
                if _rm_path(dp):
                    n += 1
                continue
            if d == "tmpfiles":
                try:
                    for child in dp.iterdir():
                        if child.name.startswith("tmp."):
                            if _rm_path(child):
                                n += 1
                except OSError:
                    pass
                continue
            if recursive:
                next_dirs.append(d)
        dirs[:] = next_dirs
        for f in files:
            fp = root_p / f
            if f.endswith(".nsf"):
                if _rm_path(fp):
                    n += 1
                continue
            if "_geom" not in f:
                continue
            if _keep_optim_file(fp):
                continue
            if f.endswith((".log", ".xyz", ".json", ".cons.inp", ".nsf")):
                if _rm_path(fp):
                    n += 1
        if not recursive:
            break
    return n


def cleanup_wavefront_geometric_scratch(
    wavefront, *, keep_incomplete_optim: bool = False
) -> None:
    """Remove leftover geomeTRIC scratch for a wavefront (not a search change)."""
    from . AseEngine import cleanup_geometric_scratch

    keep = []
    for level in getattr(wavefront, "levels", []) or []:
        for node in getattr(level, "nodes", []) or []:
            pkl = getattr(node, "node_pkl", None)
            if not pkl:
                continue
            prefix = geometric_prefix_from_node_pkl(pkl)
            keep_opt = bool(
                keep_incomplete_optim and not getattr(node, "complete", False)
            )
            cleanup_geometric_scratch(prefix, keep_optim=keep_opt)
            if keep_opt:
                keep.append(prefix)

    workdir = getattr(wavefront, "workdir", None)
    if not workdir:
        ckpt = getattr(wavefront, "checkpoint", None)
        if ckpt:
            workdir = str(Path(ckpt).resolve().parent)
    if not workdir:
        workdir = os.getcwd()
    n = sweep_geometric_scratch_dir(
        workdir, keep_optim_prefixes=keep, recursive=False
    )
    tmpfiles = Path(workdir) / "tmpfiles"
    if tmpfiles.is_dir():
        n += sweep_geometric_scratch_dir(
            tmpfiles, keep_optim_prefixes=keep, recursive=False
        )
    if n:
        print(f"[ffpopt] removed {n} leftover geomeTRIC scratch path(s) in {workdir}")
