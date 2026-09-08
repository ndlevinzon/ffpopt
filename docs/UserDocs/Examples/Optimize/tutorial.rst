.. _Optimize-tutorial:


Geometry Optimizations
======================

This tutorial will show how to perform geometry optimizations with ffpopt-Optimize.py


Learning Objectives
-------------------


Tutorial
--------


1. User Options
~~~~~~~~~~~~~~~

 .. code-block::

    usage: ffpopt-Optimize.py [-h] [-i INP] -o OUT [-n NPROC] [-m MODEL] [--mfile MFILE]
                  [--norestene] [--no-opt] [--ase-opt-tol ASE_OPT_TOL] 
		  [--psi4-num-threads PSI4_NUM_THREADS] [--psi4-memory PSI4_MEMORY]
		  [--geometric-opt] [--geometric-maxiter GEOMETRIC_MAXITER]
		  [--geometric-coordsys GEOMETRIC_COORDSYS]
		  [--geometric-converge GEOMETRIC_CONVERGE]
		  [--geometric-enforce GEOMETRIC_ENFORCE] [--geometric-ini GEOMETRIC_INI]
		  [--test] [--test-delta TEST_DELTA] [--cpu] 
    Perform a geometry optimization

    recommended options:
    ffpopt-Optimize.py --geometric-opt --geometric-ini="" --model=qdpi2 --inp=inp.json \
                       --out=out.json

    options:
      -h, --help            show this help message and exit
      -i INP, --inp INP     Input json file
      -o OUT, --out OUT     Output json file
      -n NPROC, --nproc NPROC
                            Number of optimizations to run at a time. Default 1
      --test                Finitie difference test of forces
      --test-delta TEST_DELTA
                            Finitie difference displacement (Angstroms) Default:1.e-3
      --norestene           If present, then run a single point energy after optimization
                            that excludes the restraint contributions
      -m MODEL, --model MODEL
                            Energy calculator: sander, xtb, qdpi2, dpmlp, mace, aimnet2,
			    aimnet2_wb97m (aimnet2 is an alias for aimnet2_wb97m),
			    aimnet2_b973c, aimnet2_qr, ani1x, ani2x, ani1ccx, pyscfneo,
			    theory/basis (psi4). dpmlp looks for
                            dp_test.pt in the current directory. Default: sander
      --mfile MFILE         Parameter file / network parameter file. Default: None
      --psi4-memory PSI4_MEMORY
                            The available memory. This is sent to the psi4 calculator.
			    Default: '1gb'
      --psi4-num-threads PSI4_NUM_THREADS
                            The number of threads sent to the psi4 calculator. Default: 4
      --cpu                 This will adjust the environmental variables to force machine
                            learning model evaluation on the cpus. Specifically, it exports
			    JAX_PLATFORMS='cpu' and CUDA_VISIBLE_DEVICES=-1
      --no-opt              Skip all geometry optimizations and simply perform an energy
                            evaluation
      --geometric-opt       If present, then use the BFGS optimizer in ASE rather than
                            geometric-optimize
      --ase-opt-tol ASE_OPT_TOL
                        The ASE geometry optimization tolerance. Default: 0.01 eV/A
      --geometric-maxiter GEOMETRIC_MAXITER
                        Maximum number of optimization steps. Default: 500
      --geometric-coordsys GEOMETRIC_COORDSYS
                        Coordinate system. Default: tric
      --geometric-converge GEOMETRIC_CONVERGE
                        Optimization tolerance(s). Default: 'set GAU_TIGHT'. Other named
			options include: 'set GAU_LOOSE', 'set GAU', 'set GAU_TIGHT',
			'set GAU_VERYTIGHT'
      --geometric-enforce GEOMETRIC_ENFORCE
                        Constraint enforcement tolerance. Default: 0.1
      --geometric-ini GEOMETRIC_INI
                        Path to the geometric INI file.


2. Sequential (serial) optimization of each structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 .. code-block::

    [code@user] ffpopt-Optimize.py --model=sander --geometric-opt \
                --geometric-ini="" --inp=confs.json --out=confs.sander.json


3. Threaded backfilled optimization of all structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``--nproc=5`` option will run up to 5 optimizations at the same time.

 .. code-block::

    [code@user] ffpopt-Optimize.py --model=sander --geometric-opt \
                --geometric-ini="" --inp=confs.json --out=confs.sander.json --nproc=5

		

4. MPI parallel backfilled optimization of all structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``mpirun -n 5`` option will run up to 5 optimizations at the same time, but
all 5 cores do not need to be physcailly located on the same compute node.

 .. code-block::

    [code@user] mpirun -n 5 ffpopt-Optimize.py --model=sander --geometric-opt \
                --geometric-ini="" --inp=confs.json --out=confs.sander.json

		

   
