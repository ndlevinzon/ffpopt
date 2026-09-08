"""Spawn Pool whose workers are non-daemon (may nest child pools).

Wavefront node pools are reused across sequential ``calculate()`` calls in
the same process when the worker count matches, so a multi-bond twist does
not pay spawn+import cost on every scan. Set ``FFPOPT_WF_POOL_REUSE=0`` to
create and tear down a pool per ``calculate()``.
"""

from __future__ import annotations

import atexit
import os
from typing import Any, Callable, Optional

_SpawnProcessBase = __import__("multiprocessing").get_context("spawn").Process

_NESTED_SPAWN_ENV = "FFPOPT_IN_SPAWN_WORKER"
_REUSED_POOL = {"pool": None, "key": None}
_ATEXIT_REGISTERED = False


class NonDaemonSpawnProcess(_SpawnProcessBase):
    """Spawn process that ignores daemon=True (may create nested pools)."""

    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NonDaemonSpawnContext:
    """Context wrapper so ``Pool`` uses :class:`NonDaemonSpawnProcess`."""

    def __init__(self):
        import multiprocessing as mp

        self._ctx = mp.get_context("spawn")
        self.Process = NonDaemonSpawnProcess

    def __getattr__(self, name):
        return getattr(self._ctx, name)


def in_spawn_worker() -> bool:
    """True when this process was started as an ffpopt spawn pool worker."""
    return os.environ.get(_NESTED_SPAWN_ENV, "").strip() == "1"


def _mark_spawn_worker() -> None:
    os.environ[_NESTED_SPAWN_ENV] = "1"


def _env_truthy(name: str, default: str = "0") -> bool:
    raw = os.environ.get(name, default)
    return str(raw).strip().lower() not in {"", "0", "false", "no", "off"}


def env_int(name: str, default: int = 0) -> int:
    raw = os.environ.get(name)
    if raw is None or not str(raw).strip():
        return int(default)
    try:
        return int(str(raw).strip())
    except ValueError:
        return int(default)


def pin_wavefront_worker() -> None:
    """Pool initializer: mark nested-spawn and cap math threads."""
    from . CpuThreads import pin_math_threads

    _mark_spawn_worker()
    pin_math_threads(1)


def _combined_worker_init(user_init, user_args):
    """Module-level wrapper so spawn can pickle the Pool initializer."""
    pin_wavefront_worker()
    if user_init is not None:
        user_init(*tuple(user_args or ()))


def make_nondaemon_spawn_pool(
    n_workers: int,
    *,
    initializer: Optional[Callable[..., Any]] = None,
    initargs: tuple = (),
):
    """Spawn ``Pool`` whose workers are non-daemon (may nest wavefront pools).

    ``multiprocessing.get_context(...).Pool`` is a factory method, not a class,
    so it cannot be subclassed. Pass ``multiprocessing.pool.Pool`` a context
    whose ``Process`` is :class:`NonDaemonSpawnProcess` instead.

    Workers set ``FFPOPT_IN_SPAWN_WORKER=1`` so nested bond/wavefront code can
    flatten rather than opening a third spawn level.
    """
    from multiprocessing.pool import Pool

    if initializer is None:
        init_fn = pin_wavefront_worker
        init_args: tuple = ()
    else:
        init_fn = _combined_worker_init
        init_args = (initializer, tuple(initargs))

    return Pool(
        processes=max(1, int(n_workers)),
        context=NonDaemonSpawnContext(),
        initializer=init_fn,
        initargs=init_args,
    )


def close_reused_wavefront_pool() -> None:
    """Close the process-local reused wavefront pool, if any."""
    pool = _REUSED_POOL.get("pool")
    _REUSED_POOL["pool"] = None
    _REUSED_POOL["key"] = None
    if pool is None:
        return
    try:
        pool.close()
        pool.join()
    except Exception:
        try:
            pool.terminate()
            pool.join()
        except Exception:
            pass


def _ensure_atexit() -> None:
    global _ATEXIT_REGISTERED
    if _ATEXIT_REGISTERED:
        return
    atexit.register(close_reused_wavefront_pool)
    _ATEXIT_REGISTERED = True


def acquire_wavefront_pool(
    nproc: int,
    *,
    initializer: Optional[Callable[..., Any]] = None,
    initargs: tuple = (),
    reuse_key=None,
):
    """Return ``(pool, owns_pool)`` for a wavefront ``calculate()`` drain.

    When reuse is enabled the caller must not close the pool; call
    :func:`close_reused_wavefront_pool` at workflow end (or rely on atexit).
    ``reuse_key`` should include ``nproc`` plus model / charge / parm (and any
    N-D constraint identity) once workers hold ``los`` from the initializer.
    """
    nproc = max(1, int(nproc))
    key = reuse_key if reuse_key is not None else (nproc,)
    if not _env_truthy("FFPOPT_WF_POOL_REUSE", "1"):
        return (
            make_nondaemon_spawn_pool(
                nproc, initializer=initializer, initargs=initargs
            ),
            True,
        )
    _ensure_atexit()
    if _REUSED_POOL["pool"] is not None and _REUSED_POOL["key"] == key:
        return _REUSED_POOL["pool"], False
    close_reused_wavefront_pool()
    pool = make_nondaemon_spawn_pool(
        nproc, initializer=initializer, initargs=initargs
    )
    _REUSED_POOL["pool"] = pool
    _REUSED_POOL["key"] = key
    return pool, False


def split_nproc_for_items(
    nproc: int,
    n_items: int,
    *,
    prefer_depth: bool = False,
    min_inner: int | None = None,
    flatten_nested: bool = True,
) -> tuple[int, int]:
    """Split ``nproc`` into ``(n_outer_workers, n_inner_per_worker)``.

    When ``prefer_depth`` is True, keep at least ``min_inner`` cores per outer
    worker (default from ``FFPOPT_MIN_WF_NPROC`` or 2) so HL wavefronts are not
    forced to 1-wide when many jobs share a modest allocation.

    When ``flatten_nested`` is True, never return both outer and inner greater
    than 1. Nested ``spawn`` pools (bond workers that each open a wavefront
    pool) dominate bootstrap cost inside an already-spawned worker.
    """
    nproc = max(1, int(nproc))
    n_items = max(1, int(n_items))
    if n_items == 1:
        return 1, nproc
    if not prefer_depth:
        n_outer = min(nproc, n_items)
        n_inner = max(1, nproc // n_outer)
    else:
        if min_inner is None:
            min_inner = max(1, env_int("FFPOPT_MIN_WF_NPROC", 2))
        min_inner = max(1, int(min_inner))
        max_outer = min(n_items, max(1, nproc // min_inner))
        n_outer, n_inner = 1, nproc
        best_used, best_outer = nproc, 1
        for cand_outer in range(1, max_outer + 1):
            cand_inner = nproc // cand_outer
            if cand_inner < min_inner:
                continue
            used = cand_outer * cand_inner
            if used > best_used or (used == best_used and cand_outer > best_outer):
                n_outer, n_inner = cand_outer, cand_inner
                best_used, best_outer = used, cand_outer
    if flatten_nested and n_outer > 1 and n_inner > 1:
        if prefer_depth:
            return 1, nproc
        return min(nproc, n_items), 1
    return n_outer, n_inner
