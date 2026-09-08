#!/bin/bash
set -e
set -u

# Direct geometric-optimize CLI: a new interpreter + model load per call.
# ffpopt-Optimize.py --geometric-opt uses the same converge / coordsys
# settings in-process (EngineASE + cached calculator). This directory is
# the CLI engine demo; keep geometric-optimize as the runnable command.

geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"sander","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="sander" input.xyz cons.10-9-8-6.inp
#geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"mace","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="mace" input.xyz cons.10-9-8-6.inp
#geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"qdpi2","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="qdpi2" input.xyz cons.10-9-8-6.inp
#geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"xtb","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="xtb" input.xyz cons.10-9-8-6.inp
#geometric-optimize --engine ase --ase-class="ffpopt.ase.GenCalculator" --ase-kwargs='{"mode":"HF/6-31G*","parm":"minimal_orig.parm7","crd":"minimal_orig.rst7"}' --prefix="psi4" input.xyz cons.10-9-8-6.inp

# In-process equivalent (uncomment after PrepareInput):
# python3 ../../src/python/bin/ffpopt-PrepareInput.py \
#     -p minimal_orig.parm7 -c minimal_orig.rst7 -o start.json
# python3 ../../src/python/bin/ffpopt-Optimize.py \
#     --geometric-opt -i start.json -o opt.json



