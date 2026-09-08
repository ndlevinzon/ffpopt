#!/usr/bin/env python3
"""In-process dihedral-twist example (pooled scans, reused spawn workers).

The bash CLI ``ffpopt-DihedTwistWorkflow.py`` only *emits* a driver script
and, with ``--seqscan``, uses ``ffpopt-DihedScan.py`` instead of wavefront.
This runner calls :func:`ffpopt.Workflows.run_dihed_twist_workflow` so
``nproc`` is split across independent HL/orig bond scans, workers keep the
ML calculator, and geomeTRIC runs in-process.

Must be launched under ``if __name__ == "__main__":`` (spawn multiprocessing).

Run from this directory after ``start.json`` exists (see ``run.sh``):

    python3 run.py
"""

from pathlib import Path
import sys

_LIB = Path(__file__).resolve().parents[2] / "src" / "python" / "lib"
if _LIB.is_dir():
    sys.path.insert(0, str(_LIB))

from ffpopt.Workflows import run_dihed_twist_workflow


def main():
    result = run_dihed_twist_workflow(
        inp="start.json",
        bond=["9,8", "8,6", "6,2"],
        model="qdpi2",
        geometric_opt=True,
        nproc=4,
        maxiter=2,
    )
    print(f"early_stopped_at: {result['early_stopped_at']}")
    print(f"scans: {len(result['scans'])}")
    print(f"iterations: {result['iterations']}")


if __name__ == "__main__":
    main()
