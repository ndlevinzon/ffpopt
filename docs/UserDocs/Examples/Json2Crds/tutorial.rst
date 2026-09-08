.. _json2crds-tutorial:


Extracting Coordinates
======================

This tutorial will show how to perform file conversions using ffpopt-Json2Crds.py.  This script was originally written to extract coordinates from ffpopt-style json files; however, this script can perform a number of conversions that are independent of the json file format.


Learning Objectives
-------------------


Tutorial
--------


1. Extracting all conformers from a json file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If mol.json contains multiple structures, then the following command
will output mol_s000.pdb, mol_s001.pdb, etc., where sXXX are the
structure names (the default names are sXXX where XXX is a zero-padded
integer).  If mol.json contains only 1 conformer, then the output file
would simply be mol.pdb.

.. code-block::

   [user@comp] ffpopt-Json2Crds.py --inp=mol.json --out=mol.pdb
   

Other output formats include mol2, xyz, and rst7.

   a. ".xyz" -> XYZ format
   b. ".mol2" -> MOL2 format
   c. ".rst7" -> Amber formatted restart format
   d. ".pdb" -> PDB format


2. Extracting a single conformer from a json file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can extract a single conformer by name. Suppose the name is "s000"
(the default name of the first conformer).  One can extract this conformer
with the following command. Because a single conformer is extracted,
the filename will be mol.pdb (rather than mol_s000.pdb).

.. code-block::

   [user@comp] ffpopt-Json2Crds.py --name=s000 --inp=mol.json --out=mol.pdb
   


3. Format conversion
~~~~~~~~~~~~~~~~~~~~
You do not need a json file. For example, the following converts
a pdb file to several other formats


.. code-block::

   [user@comp] ffpopt-Json2Crds.py --inp=mol.pdb --out=mol.mol2
   [user@comp] ffpopt-Json2Crds.py --inp=mol.pdb --out=mol.xyz
   [user@comp] ffpopt-Json2Crds.py --inp=mol.pdb --out=mol.rst7
   [user@comp] ffpopt-Json2Crds.py --inp=mol.pdb --out=mol.json


Note that the json file created from mol.pdb is missing information,
such as atomic charges. One would typically use ffpopt-PrepareInput.py
to create a json file, which can set constraints and restraints.


4. Inchi and Smiles strings
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following example creates an xyz file for acetic acid from its Inchi string. One can also use a smiles string.

.. code-block::

   [user@comp] ffpopt-Json2Crds.py --inp 'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)' --out mol.xyz
   [user@comp] cat mol.xyz
   8
   Properties=species:S:1:pos:R:3:initial_charges:R:1 charge=0 spin=1 pbc="F F F"
   C       -0.97938925      -0.06317125      -0.15629403       0.00000000
   C        0.46775661       0.28301764      -0.08165682       0.00000000
   O        0.85856137       1.48926781      -0.17942098       0.00000000
   O        1.46325114      -0.66733882       0.09563239       0.00000000
   H       -1.52151351       0.68733959      -0.72387465       0.00000000
   H       -1.38014251      -0.13915930       0.87249948       0.00000000
   H       -1.11880141      -1.08613953      -0.59830750       0.00000000
   H        2.21027756      -0.50381613       0.77142210       0.00000000


The next example converts an xyz file to an inchi string.


.. code-block::

   [user@comp] ffpopt-Json2Crds.py --inp mol.xyz --out tdep-inchi
   InChI=1/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/f/h3H

   [user@comp] ffpopt-Json2Crds.py --inp mol.xyz --out tindep-inchi
   InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)

   [user@comp] ffpopt-Json2Crds.py --inp mol.xyz --out tdep-smiles
   [H]OC(=O)C([H])([H])[H]

   [user@comp] ffpopt-Json2Crds.py --inp mol.xyz --out tindep-smiles
   CC(=O)O

   [user@comp] ffpopt-Json2Crds.py --inp mol.xyz --out tdep-inchitoken
   SW5DaEk9MS9DMkg0TzIvYzEtMigzKTQvaDFIMywoSCwzLDQpL2YvaDNI

   [user@comp] ffpopt-Json2Crds.py --inp SW5DaEk9MS9DMkg0TzIvYzEtMigzKTQvaDFIMywoSCwzLDQpL2YvaDNI --out tindep-smiles
   CC(=O)O

   [user@comp] ffpopt-Json2Crds.py --inp mol.xyz --out iupac
   1-hydroxy-1-oxo-ethane

   

   
The outputs are not filenames, they are special indicators of the string to output. They include tautomer-dependent and tautomer-independent variants. When creating databases, it is useful to search with tautomer-independent strings to capture all tautomers of a molecule; however, the tautomer-dependent strings are essential for rebuilding the 3d structure from the string.


   a. tdep-inchi : tautomer-dependent inchi (nonstandard inchi)
   b. tindep-inchi : tautomer-independent inchi (standard inchi)
   c. tdep-smiles : tautomer-dependent smiles
   d. tindep-smiles : tautomer-independent smiles
   e. tdep-inchikey : tautomer-dependent inchi key (output only/no back-transformation)
   f. tindep-inchikey : tautomer-independent inchi key (output only/no back-transformation)
   g. tdep-inchitoken : tautomer-dependent inchi token (ffpopt-format, URL-safe base64 tdep-inchi)
   h. tindep-inchitoken : tautomer-independent inchi token (ffpopt-format, URL-safe base64 tindep-inchi)
   i. iupac : Formal IUPAC molecule name

   
5. Charges
~~~~~~~~~~

Notice that the conversion of xyz-to-mol2 or json normally assumes the net charge of the molecule is zero, but you can set the charge by manipulating the filename to include a ":q" suffix, where q is an integer.

.. code-block::

   [user@comp] ffpopt-Json2Crds.py --inp mol.xyz:-1 --out mol.mol2
   [user@comp] grep ' MOL ' mol.mol2
   1 C           -0.9794    -0.0632    -0.1563 C             1 MOL       -0.125000
   2 C2           0.4678     0.2830    -0.0817 C             1 MOL       -0.125000
   3 O            0.8586     1.4893    -0.1794 O             1 MOL       -0.125000
   4 O2           1.4633    -0.6673     0.0956 O             1 MOL       -0.125000
   5 H           -1.5215     0.6873    -0.7239 H             1 MOL       -0.125000
   6 H2          -1.3801    -0.1392     0.8725 H             1 MOL       -0.125000
   7 H3          -1.1188    -1.0861    -0.5983 H             1 MOL       -0.125000
   8 H4           2.2103    -0.5038     0.7714 H             1 MOL       -0.125000

   
