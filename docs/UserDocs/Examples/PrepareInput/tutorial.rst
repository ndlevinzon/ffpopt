.. _PrepareInput-tutorial:


Create ffpopt Json files with ffpopt-PrepareInput.py
====================================================

This tutorial will show how to prepare a Json file used as input to
the other ffpopt scripts.


Learning Objectives
-------------------


Tutorial
--------
1. User Options
~~~~~~~~~~~~~~~

 .. code-block::

    usage: ffpopt-PrepareInput.py [-h] [--charge CHARGE] [--spin SPIN] [--parm PARM] --crd CRD [--name NAME] [--append]
           [--update] [--clear] --out OUT [--constrain CONSTRAIN] [--restrain-bond RESTRAIN_BOND]
           [--restrain-angle RESTRAIN_ANGLE] [--restrain-dihed RESTRAIN_DIHED] [--restrain-r12 RESTRAIN_R12]
	   [--restrain-puckerx RESTRAIN_PUCKERX] [--restrain-puckery RESTRAIN_PUCKERY] [--restrain-rms RESTRAIN_RMS]
           [--restrain-twist RESTRAIN_TWIST]

    Prepare na ASE XYZ file that contains charge, multiplicity, constraint, and restraint information

    options:
      -h, --help            show this help message and exit
      --charge CHARGE, -q CHARGE
                            Charge, int. Default: None (reads from parm7 or mol2 file).
      --spin SPIN, -m SPIN  Spin multiplicity, int. Default: 1
      --parm PARM, -p PARM  Amber parm7 file, only needed if crd is a rst7 file or you want to use sander in the future
      --crd CRD, -c CRD     xyz, mol2, or amber rst7 file; if --update is present, then this should be a json file
      --name NAME, -n NAME  Name of the structure. If blank, then it is sXXX, where XXX is a zero-padded number
      --append, -a          append structure to output json
      --update, -u          update an existing json file. This is used to supple a --parm value or replace the
                            constraints/restraints. If --parm is supplied and --name is not given, then all structures
			    with the same list of elements are modified. If --parm is not supplied, then only the
			    constraints/restraints are replaced (or cleared).
      --clear               if --update and --clear, then it removes all restraint and constraint definitions.
      --out OUT, -o OUT     output json file
      --constrain CONSTRAIN
                            comma-separated list of 2,3,or 4 0-based integers. The list can be appended with =value to
			    specify a value. If not given, then the input coordinates are used to assign the value.
			    This can be used multiple times.
      --restrain-bond RESTRAIN_BOND
                            str, format:'k,idx1,idx2[=value]', where idx1,idx2 are 0-based integers and value can be
			    missing
      --restrain-angle RESTRAIN_ANGLE
                            str, format:'k,idx1,idx2,idx3[=value]', where idx1,idx2,idx3 are 0-based integers and value
			    can be missing. Value should be degrees
      --restrain-dihed RESTRAIN_DIHED
                            str, format:'k,idx1,idx2,idx3,idx4[=value]', where idx1,idx2,idx3,idx4 are 0-based integers
			    and value can be missing. Value should be degrees
      --restrain-r12 RESTRAIN_R12
                            str, format:'k,idx1,idx2,idx3,idx4[=value]', where idx1,idx2,idx3,idx4 are 0-based integers
			    and value can be missing
      --restrain-puckerx RESTRAIN_PUCKERX
                            str, format:'k,idx1,idx2,idx3,idx4,idx5[=value]', where idx1,idx2,idx3,idx4,idx5 are 0-based
			    integers and value can be missing
      --restrain-puckery RESTRAIN_PUCKERY
                            str, format:'k,idx1,idx2,idx3,idx4,idx5[=value]', where idx1,idx2,idx3,idx4,idx5 are 0-based
			    integers and value can be missing
      --restrain-rms RESTRAIN_RMS
                            str, format:'[filename,]k,idx1,idx2,[...]', where idx1,idx2,[...] are 0-based integers. If
			    filename is missing, then the input coordinates are adopted. If filename is present, it
			    must be a XYZ file. The indexes reflect those atoms that contribute as nonzero mass-weighted
			    RMS contribution.
      --restrain-twist RESTRAIN_TWIST
                            str, format:'[filename,]k,idx1,idx2', where idx1,idx2 are 0-based integers. If filename is
			    missing, then the input coordinates are adopted. If filename is present, it must be a XYZ
			    file. The indexes reflect the bond that is being twisted, which will create separate rms
			    restraints for each half of the bond.


2. Create a json from various input formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prepare from an Amber parm7/rst7 so one can run sander MM geometry optimizations.  The json file format has a "parm" value that is the filename of the parm7 file used to initialize pysander.

.. code-block::

   [user@comp] ffpopt-PrepareInput.py --parm mol.parm7 --crd mol.rst7 --out mol.json


One can also generate the json from a mol2 file, but one won't be able to run
sander optimizations unless you fill-in the missing "parm" value.
   
.. code-block::

   [user@comp] ffpopt-PrepareInput.py --crd mol.mol2 --out mol.json


Other input formats are possible. You can use the --charge option to specify a net charge.  The atom types will simply be the element symbols. The charges will be netcharge/numatoms.  If all you have is an inchi/smiles/xyz etc, then you'd be better off using ffpopt-antechamber.py to assign am1-bcc charges and gaff atom types, for example.

.. code-block::

   [user@comp] ffpopt-PrepareInput.py --charge=0 --crd mol.xyz --out mol.json


3. Adding Restraints
~~~~~~~~~~~~~~~~~~~~

a. -\-restrain-bond="k,idx1,idx2[=value]"
   
   i. Restrain bond between atoms 0 and 1 using 50 eV/A^2 force constant
   with a harmonic potential centered about the bond length measured from
   the input coordinates.

   ``--restrain-bond="50.0,0,1"``

   ii. Restrain bond between atoms 0 and 1 using 50 eV/A^2 force constant
   with a harmonic potential centered about 2.1 A.

   ``--restrain-bond="50.0,0,1=2.1"``
   
b. -\-restrain-angle="k,idx1,idx2,idx3[=value]"
   
   i. Restrain angle between atoms 0-1-2 using 50 eV/A^2 force constant
   with a harmonic potential centered about the angle measured from the
   input coordinates.

   ``--restrain-angle="50.0,0,1,2"``

   ii. Restrain angle between atoms 0-1-2 using 50 eV/A^2 force constant
   with a harmonic potential centered about 30 degrees.
       
   ``--restrain-angle="50.0,0,1=30."``

c. -\-restrain-dihedral="k,idx1,idx2,idx3,idx4[=value]"
   
   i. Restrain angle between atoms 0-1-2-4 using 50 eV/A^2 force constant
   with a harmonic potential centered about the angle measured from the
   input coordinates.

   ``--restrain-dihed="50.0,0,1,2"``

   ii. Restrain angle between atoms 0-1-2-4 using 50 eV/A^2 force constant
   with a harmonic potential centered about 30 degrees.
       
   ``--restrain-dihed="50.0,0,1=30."``

d. -\-restrain-r12="k,idx1,idx2,idx3,idx4[=value]"

   i. Restrain r12 =
   R(idx1,idx2) - R(idx3,idx4)

e. -\-restrain-puckerx="k,idx1,idx2,idx3,idx4,idx5[=value]"

   i. Indexes 1-5 coorespond to C1', C2', C3', C4', O4'.
   This is the X-coordinate from
   `doi.org/10.1021/ct401013s <https://doi.org/10.1021/ct401013s>`_.

f. -\-restrain-puckery="k,idx1,idx2,idx3,idx4,idx5[=value]"

   i. Indexes 1-5 coorespond to C1', C2', C3', C4', O4'.
   This is the Y-coordinate from
   `doi.org/10.1021/ct401013s <https://doi.org/10.1021/ct401013s>`_.

g. -\-restrain-rms="[filename,]k,idx1,idx2,[...]"

   i. If filename is missing, then the input coordinates are adopted.
   If filename is present, it must be a XYZ file.
   The indexes reflect those atoms that contribute as nonzero
   mass-weighted RMS contribution.

h. -\-restrain-twist='[filename,]k,idx1,idx2'

   i. If filename is missing, then the input coordinates are adopted.
   If filename is present, it must be a XYZ file.
   The indexes reflect the bond that is being twisted, which will
   create separate rms restraints for each half of the bond.

      
4. Adding Constraints
~~~~~~~~~~~~~~~~~~~~~
Constraints only work well if you have a single constraint.
If you need to maintain several conditions, you should prefer restraints.

a. Constrain a bond.
   ``--constrain="idx1,idx2"``

b. Constrain an angle.
   ``--constrain="idx1,idx2,idx3"``

c. Constrain an dihedral.
   ``--constrain="idx1,idx2,idx3,idx4"``

   
5. Clear all constraints and restraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

   [user@comp] ffpopt-PrepareInput.py --update --clear --crd mol.json --out cleared.json

   
6. Clear constraints and restraints for a specific molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

   [user@comp] ffpopt-PrepareInput.py --name="s000" --update --clear --crd mol.json --out cleared.json

   
7. Append a conformer to an existing json file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

   [user@comp] ffpopt-PrepareInput.py --append --crd other.mol2 --out mol.json


   
   


