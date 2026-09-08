#!/bin/bash
set -e
set -u


# ``--geometric-opt`` now runs geomeTRIC in-process against a cached ASE
# calculator (set FFPOPT_GEOMETRIC_SUBPROCESS=1 for the old CLI path).
# ``--parallel`` spawn-pools the forward and reverse scans (nondaemon workers).

python3 ../../src/python/bin/ffpopt-PrepareInput.py \
	--parm minimal_orig.parm7 \
	--crd minimal_orig.rst7 \
	--out start.json

python3 ../../src/python/bin/ffpopt-DihedScan.py \
	--inp start.json \
	--out scan.json \
	--dihed="10,9,8,6" \
	--delta=15 --parallel --geometric-opt

python3 ../../src/python/bin/ffpopt-Json2Crds.py \
        --inp scan.json \
        --out scan.xyz

