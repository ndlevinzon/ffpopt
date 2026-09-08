#!/bin/bash
set -e
set -u

tdir=../../src/python/bin

python3 ${tdir}/ffpopt-PrepareInput.py \
	-p minimal_orig.parm7 -c minimal_orig.rst7 \
	-o start.json

python3 ${tdir}/ffpopt-Optimize.py \
	-i start.json -o oscan.json
# ASE BFGS by default (persistent calculator cache). Add --geometric-opt
# for in-process geomeTRIC; add -n N for ParallelGeomOpt spawn workers.


python3 ${tdir}/ffpopt-Json2Crds.py \
	-i oscan.json -o oscan.xyz




