.. _antechamber-tutorial:


Create mol2 Files with ffpopt-antechamber.py
============================================

This tutorial will show how to prepare a mol2 file prepared using antechamber via the ffpopt-antechamber.py script.
This script allows multiple input formats (json, mol2, xyz, pdb, smiles, inchi).
It provides options too perform simply conformational searches and conformational averaging of the charges.
It also provides an easy-to-use machanism to prevent behind-the-scenes am1 optimizations with sqm.



Learning Objectives
-------------------


Tutorial
--------


1. User Options
~~~~~~~~~~~~~~~

 .. code-block::

    usage: ffpopt-antechamber.py [-h] -i INP -o OUT [-c {bcc,am1bcc,abcg2}]
           [-at {gaff,gaff2,amber,bcc,abcg2,sybyl}] [-rn RESIDUE_NAME] [--verbose] [--clean]
	   [--nconf NCONF] [--nkeep NKEEP] [--no-opt] [--confdep]

    ffpopt parameterization wrapper mapping traditional antechamber CLI flags to RunCreateAmberMol2.

    options:
      -h, --help            show this help message and exit
      -i INP, --inp INP, --input-file INP
                            Path to the molecular structure input target (e.g., .xyz, .pdb,
			    .mol2, .json, or inchi/smiles). If the input is an xyz or pdb,
			    then the name can be suffixed with a charge. For example:
			    --inp=foo.xyz:-1 means the net charge of the molecule is -1
			    --inp=foo.pdb:1 means the net charge of the molecule is +1.
			    The default is to assume the molecule is neutral.
      -o OUT, --out OUT, --output-file OUT
                            Output .mol2 or .json
      -c {bcc,am1bcc,abcg2}, --charge-method {bcc,am1bcc,abcg2}
                            Antechamber charge evaluation routing scheme choice. Complete
			    valid choices: bcc, am1bcc, abcg2. Default: am1bcc.
      -at {gaff,gaff2,amber,bcc,abcg2,sybyl},
      --atom-type {gaff,gaff2,amber,bcc,abcg2,sybyl}
                            Force field atom typing family classification matrix rules. Complete
			    valid choices: gaff, gaff2, amber, bcc, abcg2, sybyl. Default: gaff.
      -rn RESIDUE_NAME, --resname RESIDUE_NAME
                            Residue designator code naming convention token to label the system
			    outputs. Default: MOL.
      -v, --verbose         Turn on verbose printing.
      --clean               Purge intermediate runtime file wrappers and scratch environments
                            when finished.
      --nconf NCONF         Number of initial candidate configuration geometries to generate
                            inside the conformer pool search loop. Set to 0 to skip entirely.
			    Default: 0
      --nkeep NKEEP         Maximum structural tracking storage threshold count used to clamp
                            the kept lowest-energy ensemble pool members. Default: 100
      --no-opt              If present, don't perform a am1 optimization when calculating charges
      --confdep             If present, save conformer-dependent charges to {base}.confdepqs.{ext}
                            where base and ext are the base name and extension of --out



2. GAFF example
~~~~~~~~~~~~~~~
This example runs a pdb structure through antechamber to assign gaff atom types and am1-bcc charges. The name of the residue will be "MOL".  The output is mol.json.

.. code-block::

   [user@comp] ffpopt-antechamber.py --inp=inp.pdb --out mol.json -c am1bcc -at gaff -rn MOL --clean

   
3. GAFF2 example from inchi string
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example constructs a 3d structure from an inchi string, performs a conformer search, exacts the lowest energy conformation, and runs the structure through antechamber to assign gaff2 atom types and abcg2 charges. The name of the residue will be "MOL".  The output is mol.json.

.. code-block::

   [user@comp] ffpopt-antechamber.py --inp='InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3' \
   --out mol.json -c abcg2 -at gaff2 -rn MOL --nconf=100 --nkeep=1 --clean


4. GAFF2 example with conformational averaging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example is similar to the previous, except (up to) 5 of the lowest-energy conformers are extracted frin the search. Antechamber is told not to perform a am1 optimization. Instead a am1 single point calculation is performed with each structure, and mol.json contains the conformer-averaged charges.

.. code-block::

   [user@comp] ffpopt-antechamber.py --inp='InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3' \
   --out mol.json -c abcg2 -at gaff2 -rn MOL --nconf=100 --nkeep=5 --no-opt --clean


5. Conformer-dependent abcg2 charges
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example is similar to the previous, the ``--confdep`` option is used.
This causes a special file mol.confdepqs.json (the output extension is changed from .json to .confdepqs.json).  The structures within this special file are the individual charges obtained from each conformer in isolation.

.. code-block::

   [user@comp] ffpopt-antechamber.py --inp='InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3' \
   --out mol.json -c abcg2 -at gaff2 -rn MOL --nconf=100 --nkeep=5 --no-opt --confdep --clean



   
   
