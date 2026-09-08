.. _respfit-tutorial:


Restrained Electrostatic Potential Charge Fitting
=================================================

This tutorial will show how to perform RESP calculations to derive MM charges.


Learning Objectives
-------------------


Tutorial
--------

1. User Options
~~~~~~~~~~~~~~~

 .. code-block::
    
    usage: ffpopt-RespFit.py [-h] --out OUT --inp INP [--respf] [--espaloma] [--hilfiker] [--nofit]
           [--program PROGRAM] [--resp-a RESP_A] [--resp-b RESP_B] [--density DENSITY]
	   [--digits DIGITS] [--scosmo SCOSMO] [--ext-scale EXT_SCALE]
           [--ext-density EXT_DENSITY] [--group GROUP] [--freeze FREEZE] [--update-only-grouped-atoms]
	   [--update-only-ungrouped-atoms] [--confdep] [-m MODEL] [--mfile MFILE]
	   [--psi4-memory PSI4_MEMORY] [--psi4-num-threads PSI4_NUM_THREADS] [--cpu]
           confs [confs ...]

    Read or create a structure and search for conformations

    positional arguments:
      confs                 1-or-more conformers. Either json, xyz, or mol2 files

    options:
      -h, --help            show this help message and exit
      --out OUT             Output json or mol2 file.
      --inp INP             Input json or mol2 file. The output will be the same as this file but
                            with different charges. If this is a json file, then only the first
			    structure is examined.
      --respf               If present, perform resp fit using the 'resp' program (expected to be
                            in path)
      --espaloma            If present, then use conformer-independent espaloma charges from
                            https://doi.org/10.1021/acs.jpca.4c01287
      --hilfiker            If present, then use conformer-dependent hilfiker charges from
                            https://doi.org/10.48550/arXiv.2512.13579
      --nofit               If present, then charge fitting is not performed; instead the input
                            charges are rounded to the desired digits and the forced to sum to the
			    nearest integer
      --program PROGRAM     Ab initio executable. Default: psi4. This could also be gaussian;
                            e.g., --program=g16
      --resp-a RESP_A       Hyperbolic penalty prefactor. Default: 0.001. The penalty is
                            pen= a * sum_i ( sqrt( qi**2 + b**2 ) - b ), where i loops over all
			    heavy atoms.
      --resp-b RESP_B       Hyperbolic penalty width. Default: 0.1. The penalty is
                            pen= a * sum_i ( sqrt( qi**2 + b**2 ) - b ), where i loops over all
			    heavy atoms.
      --density DENSITY     The density of surface points. Default: 6 pts/Ang**2.
      --digits DIGITS       Round charges to this number of digits. Default: 4
      --scosmo SCOSMO       Calculate fits in gas and cosmo environments and take a weighted
                            average of the 2 charge vectors. The SCOSMO perturbation is a response
			    to the MM charges fit to the gas phase ESP. Default: 0.0, which
			    returns only the gas phase charges.
      --ext-scale EXT_SCALE
                            UFF radius scale factor used to generate the external potential
			    surface. Default: 1.1
      --ext-density EXT_DENSITY
                            Density of external potential points. Default: 2 pts/Ang**2
      --group GROUP         Amber-style atom selection string indicating a group of atoms whose
                            charge-sum should not change. This can be used multiple times. This
			    could potentially fail if an atom exists in multiple groups.
      --freeze FREEZE       Amber-style atom selection string indicating a group of atoms whose
                            charge should not change. This is the same as creating a --group for
			    each atom in the selection
      --update-only-grouped-atoms
                            If present, then the output file only update the atomic charges for
			    those atoms present in any group (or freeze). This can only be used
			    if atoms are not part of multiple groups.
      --update-only-ungrouped-atoms
                            If present, then the output file only update the atomic charges for
			    those atoms not present in any group (or freeze). This can only be
			    used if atoms are not part of multiple groups.
      --confdep             If present, write a file of conformer-dependent charges. The name
                            of the file is the same as 'out' but the {ext} extension is replaced
			    with .confdepqs.{ext}
      -m MODEL, --model MODEL
                            Energy calculator: sander, xtb, qdpi2, dpmlp, mace, aimnet2,
			    aimnet2_wb97m(aimnet2 is an alias for aimnet2_wb97m), aimnet2_b973c,
			    aimnet2_qr, ani1x, ani2x, ani1ccx, pyscfneo, theory/basis (psi4).
			    dpmlp looks for dp_test.pt in the current directory. Default: sander
      --mfile MFILE         Parameter file / network parameter file. Default: None
      --psi4-memory PSI4_MEMORY
                            The available memory. This is sent to the psi4 calculator.
			    Default: '1gb'
      --psi4-num-threads PSI4_NUM_THREADS
                            The number of threads sent to the psi4 calculator. Default: 4
      --cpu                 This will adjust the environmental variables to force machine
                            learning model evaluation on the cpus. Specifically, it exports
			    JAX_PLATFORMS='cpu' and CUDA_VISIBLE_DEVICES=-1



2. Create the molecule using antechamber 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following reads an inchi string,
performs a conformer search,
clusters the results,
and saves up to 5 of the lowest energy conformations.
The conformations are internally optimized by RDKit using the mmff94 force field.
It then runs each conformer through antechamber.
The conformations are written as mol2 files so you can view them.

.. code-block::

   [user@comp] ffpopt-antechamber.py --inp='InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3' \
               --out mol.json -c abcg2 -at gaff2 -rn MOL --nconf=100 --nkeep=5 \
	       --clean --verbose
   [user@comp] ffpopt-Json2Crds.py --inp mol.json --out tmp.mol2
   [user@comp] ls tmp*.mol2
   >>> tmp_s000.mol2  tmp_s001.mol2  tmp_s002.mol2

   
The clustering algorithm found 3 conformations.
The script runs antechamber on each of the conformers to obtain charges for each configuration.
The above code is similar to the commands shown below, except we keep up to 5 conformers
and run antechamber on each of the conformers (instead of keeping 1 conformer and
running antechamber 1 time).
  

.. code-block::

   [user@com] ffpopt-ConfSearch.py --out=mol.json --nkeep=1 'InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3'
   [user@com] ffpopt-Json2Crds.py --inp=mol.json --out=mol.pdb
   [user@com] antechamber -i mol.pdb -fi pdb -o am1bcc.mol2 -fo mol2 -c abcg2 -nc 0 -rn MOL \
              -at gaff2 -du y -an y -pf y -seq n


In order to do a RESP fit, we need to run the molecule through antechamber first because
the ffpopt-RespFit.py script does not try to automatically determine atoms that should have
the same charge (for example, all 3 hydrogens on a methyl group should have the same charge).
Instead, it follows the pattern that it reads from the input molecule; if the input molecule
has atoms with the same atomic number and same charge, then it will enforce this condition
on the output RESP-fitted charges.
If you already have a PDB file of your molecule, you should still run it through antechamber
to create am1-bcc (or abcg2) charges for the sole purpose of identifying the charge
symmetries within the molecule.


3. Perform the ab initio calculations and the RESP fit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ffpopt-RespFit.py script can use several ab initio backends to perform the necessary SCF calculations.
The default is the use psi4, but this can be changed with the ``--program`` option.
The script will first check to see if the ab initio input (*_*.inp) and output (*_*.log) files are present.
If the files are present, then it reads the output without running the ab initio calculation to save time.
It does this without checking if their contents are consistent with your command-line arguments; therefore,
you should delete (or move) these files before rerunning the RESP fit with different options.

.. code-block::

   [user@comp] ffpopt-RespFit.py --model="hf/6-31g*" --inp=mol.json --out=resp.json mol.json
   [user@comp] rm *_*.inp *_*.log
   [user@comp] rm mol.json

One can extract the mol2 file containing the multi-conformer RESP fit charges using ffpopt-Json2Crds.py.
The usage shown here extracts the first conformer (the lowest energy conformer).

.. code-block::

   [user@comp] ffpopt-Json2Crds.py --name=s000 --inp=resp.json --out=resp.mol2
   
   
4. Find out more
~~~~~~~~~~~~~~~~
This short example does not highlight all of the command line options.
For the most up-to-date list of options, run ``ffpopt-RespFit.py --help``.
Some of the most used options are listed below.

   a. ``--program="command"``. This sets the command used to run the ab initio softwarem, where the "command" could be "g16", "mpirun -n 4 quick.MPI", "quick.cuda", "psi4", etc.  It supports Gaussian, Psi4, and quick (only for gas phase fits).
   b. ``--respf``. If present, then perform the fit using Kollman's resp.f program rather than ffpopt's internal optimization procedure based on the description in 10.1021/j100142a004
   c. ``--group="mask"``. This option can be used multiple times to enforce constraints on the charges. The sum of fitted charges within a group will be preserved from the input. The mask is an Amber mask; e.g., "@C9,H91,H92,H93".
   d. ``--scosmo=float``. If present, then perform the fit in the gas phase and a second fit in a "COSMO-like" environment and take a linear combination of the two charge arrays.  The default is 0.0, which corresponds to a gas phase electrostatic potential.
   e. If you are using the pytorch version of ffpopt, you can use the ``--espaloma`` or ``--hilfiker`` models. The espamola model is conformer independent.
   f. The ``--nofit`` options mere rounds the charges to the desired digits and the forces the sum to round to the nearest integer.

      
