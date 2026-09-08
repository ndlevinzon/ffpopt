#!/bin/bash
set -e
set -u

# Prefer this checkout's ffpopt so in-process wavefront / cached ASE apply
# even if another copy is on PYTHONPATH.
_FFPOPT_LIB="$(cd "$(dirname "$0")/../.." && pwd)/src/python/lib"
export PYTHONPATH="${_FFPOPT_LIB}${PYTHONPATH:+:$PYTHONPATH}"

if [ ! -e start.json ]; then
    python3 ../../src/python/bin/ffpopt-PrepareInput.py \
	    --parm minimal_orig.parm7 --crd minimal_orig.rst7 \
	    --out start.json
fi

# In-process twist (drop-in for the performance ports): pooled HL+orig
# wavefronts, reused spawn workers, in-process geomeTRIC. ``nproc`` is split
# across the three bonds. Requires the ``if __name__ == "__main__"`` guard
# in run.py.
python3 run.py

# Legacy: emit a bash driver (one new process per bond; ``--seqscan`` uses
# DihedScan, not wavefront). Uncomment to compare against the old path.
# ../../src/python/bin/ffpopt-DihedTwistWorkflow.py \
#     --geometric-opt --nproc=4 \
#     --inp start.json \
#     --bond=9,8 --bond=8,6 --bond=6,2 \
#     --model=qdpi2 --seqscan > scan.sh
# bash ./scan.sh
