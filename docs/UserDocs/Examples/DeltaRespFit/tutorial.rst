Delta RESP Charge Fitting for Mutated Standard Residues
=======================================================

This tutorial will show how to derive charges for a mutated
standard residue. The specific example creates a protonated
adenine nucleobase at the N3 position.


The general idea is to perform a RESP fit of the native molecule,
a RESP fit of the mutated molecule, and add the difference
between these 2 charges sets to the native residue's
standard Amber charges.


The steps involved in this procedure are nearly identical
to the description in the ffpopt-RespFit.py tutorial
:ref:`respfit-tutorial`.


Learning Objectives
-------------------


Tutorial
--------

1. User Options
~~~~~~~~~~~~~~~

 .. code-block::

    usage: ffpopt-DeltaRespFit.py [-h] --native NATIVE --modified MODIFIED --out OUT [--respf]
           [--espaloma] [--hilfiker] [--program PROGRAM] [--resp-a RESP_A] [--resp-b RESP_B]
	   [--density DENSITY] [--digits DIGITS] [--scosmo SCOSMO]
           [--ext-scale EXT_SCALE] [--ext-density EXT_DENSITY] [--native-cap NATIVE_CAP]
	   [--modified-cap MODIFIED_CAP] [--mcss-allow-mismatch] [--dont-modify-cap-charges]
	   [--mcss-map MCSS_MAP] [-m MODEL]
           [--mfile MFILE] [--psi4-memory PSI4_MEMORY] [--psi4-num-threads PSI4_NUM_THREADS]
	   [--cpu]

    Perform a geometry optimization

    options:
    -h, --help            show this help message and exit
    --native NATIVE       The json (or mol2) file describing the native system. If it is a json,
                          then all structures within the json are assumed to be conformations
			  used to perform a multi-conformer fit.
    --modified MODIFIED   The json (or mol2) file describing the native system. If it is a json,
                          then all structures within the json are assumed to be conformations
			  used to perform a multi-conformer fit.
    --out OUT             Output json or mol2 file.
    --respf               If present, perform resp fit using the 'resp' program (expected to be
                          in path)
    --espaloma            If present, then use conformer-independent espaloma charges from
                          https://doi.org/10.1021/acs.jpca.4c01287
    --hilfiker            If present, then use conformer-independent hilfiker charges from
                          https://doi.org/10.48550/arXiv.2512.13579
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
                          average of the 2 charge vectors. The SCOSMO perturbation is a
			  response to the MM charges fit to the gas phase ESP. Default: 0.0,
			  which returns only the gas phase charges.
    --ext-scale EXT_SCALE
                          UFF radius scale factor used to generate the external potential
			  surface. Default: 1.1
    --ext-density EXT_DENSITY
                          Density of external potential points. Default: 2 pts/Ang**2
    --native-cap NATIVE_CAP
                          Amber-style atom selection string indicating a group of atoms whose
			  charge-sum should not change. This can be used multiple times. This
			  could potentially fail if an atom exists in multiple groups.
    --modified-cap MODIFIED_CAP
                          Amber-style atom selection string indicating a group of atoms whose
			  charge-sum should not change. This can be used multiple times. This
			  could potentially fail if an atom exists in multiple groups.
    --mcss-allow-mismatch
                          Default does not allow MCSS mapping between different elements. If
			  present, then a carbon can map to a flourine, for example
    --dont-modify-cap-charges
                          If present, then the atomic charges in the caps are the exact same
			  as the input. The default is to replace the cap charges with the
			  result of the constrained RESP fit. In both cases, the net charge
			  of each cap is preserved.
    --mcss-map MCSS_MAP   File containing lines like: AT1 => AT2 where AT1 is an atom name in
                          native and AT2 is an atom name in 2. If present, then this mapping
			  is used instead of MCSS
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



2. Create mol2 files for the native and modified residues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example starts with mol2 files for 
:download:`adenine (A.mol2) <A.mol2>`
and a
:download:`N3-protonated adenine (A3.mol2) <A3.mol2>`.
These files were prepared with tleap.
In both files, the nucleotide sugar has been replaced
by a methyl group (atom selection "@C9,H91,H92,H93")
connected to the N7 position.
This methyl group shall be referred to as the "cap".
A cap doesn't have to be a methyl group; it is merely
a small functional group used to take the place of
a much larger subunit.
The net charge of the methyl group (0.1241) was chosen
to neutralize the native nucleobase.
The remaining charges within A3.mol2 are the same as A.mol2;
our goal is to figure out what those charges should
be.  The exception is the charge of the H3 proton,
which was added to the A3 residue at the N3 position.
The H3 atom has been given a charge a full 1+ charge,
such that the protonated N3 residue has a net charge of 1+.



3. Create 1-or-more conformers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RESP fitting is more robust if multiple conformers are considered during the fit.  We can use the ffpopt-ConfSearch script to search for conformers and output them to confs.json.  More information can be found in the conformer search tutorial:  :ref:`confsearch-tutorial`.
For Delta RESP fits, we look for conformers of both molecules.

.. code-block::

   [user@comp] ffpopt-ConfSearch.py --out=Aconfs.json A.mol2
   [user@comp] ffpopt-ConfSearch.py --out=A3confs.json A3.mol2


4. Perform the ab initio calculations and the RESP fit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

   [user@comp] ffpopt-DeltaRespFit.py \
       --model="hf/6-31g*" \
       --native=Aconfs.json \
       --modified=A3confs.json \
       --out=A3resp.json \
       --native-cap="@C9,H91,H92,H93" \
       --modified-cap="@C9,H91,H92,H93"
   [user@comp] rm *_*.inp *_*.log
   [user@comp] rm Aconfs.json A3confs.json


The script produces a lot of output to the screen; however, the important
thing is the output file A3resp.json.
The screen will display several items, described below.


The following lists the caps for native (Cap1) and the modified (Cap2)
residues and their charge sum (q).

.. code-block::

   Cap definition consistency between native and modified residues
   Cap1:      @C9,H91,H92,H93 q: 0.124100   Cap 2:      @C9,H91,H92,H93 q: 0.124100



The following lists the "softcore" atoms from the native (sc1) and modified (sc2) residues. By "softcore", it is meant that the sc1 atoms are not present in the modified residue, and the sc2 atoms are not present in the native residue

.. code-block::
   
   MCSS Softcore selections
   sc1: []
   sc2: ['H3']


The following associates equivalent (0-based) atom indexes between the native and modified residues.


.. code-block::
   
   MCSS Common-core mapping (non-capped atoms)
   12 =>   12
   13 =>   13
    4 =>    4
    5 =>    5
    9 =>    9
   10 =>   10
   11 =>   11
    6 =>    6
    7 =>    7
    8 =>    8
    3 =>    3
    1 =>    1
    0 =>    0
    2 =>    2


The following output verifies that the sum of cap charges are the same before and after the RESP fits.


.. code-block::
   
   Conservation of cap charges (native)
   orig sum=  0.12410000  fit sum=  0.12410000

   Conservation of cap charges (modified)
   orig sum=  0.12410000  fit sum=  0.12410000



The following shows how the charges of the modified residue are calculated.
The first column is the atom name from the modified residue.
o1 is the original charge from the native residue before fitting.
f1 is the RESP fit charge from the native residue.
f2 is the RESP fit charge from the modified residue.
qnew is o1+(f2-f1).
The final charge is qnew after rounding to the desired number of decimal places (typically 4).
Below the atom charges, it provides the sum of charges in the rows above, the unseen cap
charges (which remain unchanged because one normally throws the cap away after fitting),
the charge contribution from the softcore atoms not listed in the above rows, and
the net residue charge (sum+cap+sc).

.. code-block::
   
   Charge changes of non-capped atoms
   N9 o1= -0.025100 f1= -0.288879 f2= -0.160021 qnew=  0.103758 final=  0.103800
   C8 o1=  0.200600 f1=  0.222182 f2=  0.241453 qnew=  0.219871 final=  0.219900
   H8 o1=  0.155300 f1=  0.169123 f2=  0.218968 qnew=  0.205145 final=  0.205100
   N7 o1= -0.607300 f1= -0.537340 f2= -0.482060 qnew= -0.552020 final= -0.552000
   C5 o1=  0.051500 f1=  0.016836 f2= -0.139713 qnew= -0.105049 final= -0.105000
   C6 o1=  0.700900 f1=  0.701393 f2=  0.909247 qnew=  0.908754 final=  0.908800
   N6 o1= -0.901900 f1= -0.889897 f2= -0.762345 qnew= -0.774348 final= -0.774300
   H61 o1=  0.411500 f1=  0.373926 f2=  0.382056 qnew=  0.419630 final=  0.419600
   H62 o1=  0.411500 f1=  0.373926 f2=  0.382056 qnew=  0.419630 final=  0.419600
   N1 o1= -0.761500 f1= -0.670292 f2= -0.696095 qnew= -0.787302 final= -0.787300
   C2 o1=  0.587500 f1=  0.415953 f2=  0.631561 qnew=  0.803108 final=  0.803100
   H2 o1=  0.047300 f1=  0.116399 f2=  0.149925 qnew=  0.080826 final=  0.080800
   N3 o1= -0.699700 f1= -0.624746 f2= -0.705621 qnew= -0.780575 final= -0.780600
   C4 o1=  0.305300 f1=  0.497318 f2=  0.509560 qnew=  0.317541 final=  0.317500
   H3 o1=  0.000000 f1=  0.000000 f2=  0.396930 qnew=  0.396930 final=  0.396900
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   sum o1= -0.124100 f1= -0.124100 f2=  0.875900 qnew=  0.875900 final=  0.875900
   cap o1=  0.124100 f1=  0.124100 f2=  0.124100 qnew=  0.124100 final=  0.124100
   sc o1=  0.000000 f1=  0.000000 f2=  0.000000 qnew=  0.000000 final=  0.000000
   net o1=  0.000000 f1= -0.000000 f2=  1.000000 qnew=  1.000000 final=  1.000000


   
5. Extract the mol2 file
~~~~~~~~~~~~~~~~~~~~~~~~

One can extract the mol2 file containing the multi-conformer RESP fit charges using ffpopt-Json2Crds.py/

.. code-block::

   [user@comp] ffpopt-Json2Crds.py --inp=A3resp.json --out=A3resp.mol2
   
   
6. Find out more
~~~~~~~~~~~~~~~~
This short example does not highlight all of the command line options.
For the most up-to-date list of options, run ``ffpopt-DeltaRespFit.py --help``.
Some of the most used options are listed below.

   a. ``--program="command"``. This sets the command used to run the ab initio softwarem, where the "command" could be "g16", "mpirun -n 4 quick.MPI", "quick.cuda", "psi4", etc.  It supports Gaussian, Psi4, and quick (only for gas phase fits).
   b. ``--respf``. If present, then perform the fit using Kollman's resp.f program rather than ffpopt's internal optimization procedure based on the description in 10.1021/j100142a004
   c. ``--scosmo=float``. If present, then perform the fit in the gas phase and a second fit in a "COSMO-like" environment and take a linear combination of the two charge arrays.  The default is 0.0, which corresponds to a gas phase electrostatic potential.
