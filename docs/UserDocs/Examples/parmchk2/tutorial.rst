.. _parmchk2-tutorial:


Create frcmod Files with ffpopt-parmchk2.py
===========================================

This tutorial will show how to create a frcmod file using ffpopt-parmchk2.py.
This script is capable of reading malformed parameter files, such as modxna.
Furthermore, the output frcmod file will contain all of the MASS and NONB
entries needed for parmed to load the file correctly.
Finally, the script uses the gaff or gaff2 force field to fill-in missing parameters
labeled ATTN.

This utility reads a molecular structure or file (such as a Mol2 or SMILES string),
extracts structural connectivity metrics, and executes parmchk2. It outputs a
Force Field Modification (.frcmod) file containing missing parameter definitions 
reconstructed with reasonable parameter suggestions.



Learning Objectives
-------------------


Background
----------

1. Force field type
~~~~~~~~~~~~~~~~~~~

One prepares a molecule with atom types from a "base force field",
typically either gaff, gaff2, parm99, parm10, or lipid14.
We shall call this the "force field type".
**The force field type is the "-s" command line argument when using parmchk2
(-s or -\-fftype when using ffpopt-parmchk2.py).**


2. Standard force field
~~~~~~~~~~~~~~~~~~~~~~~

Given a set of atom types, one loads a set of standard set of
force field parameters that supplement or override the base force field.
Examples include:

  1. Protein: ff99SB, ff14SB, ff03
  2. DNA: ol15, bsc1
  3. RNA: ol3, yil

These are often combined, for example "ff14SB+ol15+ol3".
**The combination of standard force fields is the "-frc" command line
argument when using parmchk2 (-frc or -\-ffstd when using ffpopt-parmchk2.py).**


3. User force field
~~~~~~~~~~~~~~~~~~~

The user may have their own parameter database. An explicit example
of this is `modxna.frcmod <https://modxna.chpc.utah.edu/>`_.
The user force field is often not used except is special situations.
Instead of having a special name, it is simply an input frcmod file.
**The user force field is the "-afrc" command line argument when using
parmchk2 (-afrc or -\-ffuser when using ffpopt-parmchk2.py).**



4. Missing-parameter force field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When running parmchk2, it may have difficulty assigning some of the
force field parameters. In these cases, it will print "ATTN" next
to the parameter within the output frcmod file.
One can use the gaff or gaff2 force fields to look-up values for
these parameters.
We call this the "missing parameter force field".
**The missing parameter force field is the "-\-ffmiss" command line argument when using
ffpopt-parmchk2.py.**


5. Atom scoring of new  atom types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parmchk2 will search for parameters in the base,
standard, user force field definitions.
If it doesn't exist, it will try to find
similar force field parameters using an "atom type score".
The scores for atoms within the base and standard force fields
are located within ${AMBERHOME}/dat/antechamber/PARMCHK.DAT.
If you create new atom types, you will either need to modify
that global score file or create a supplementary file.
**The supplemenary score file is the "-c" option when using
parmchk2 (-c or -\-score-file when using ffpopt-parmchk2.py).**


To fully add a custom atom type, you must append an independent multi-line block using this exact column layout and capitalization:

 .. code-block::

    PARM    [atomtype]    [improper_flag]    [group_id]    [mass]    [equivalent_flag]    [atomic_num]
    EQUA    [gaff_equivalent]
    CORR    [related_type_A]   [penalty matrix weights...]

    
* PARM Line Fields (Whitespace-separated columns):
    * Column 1: PARM (Literal index string)
    * Column 2: ``[atomtype]`` -- Your custom 1-to-2 character atom type identifier (e.g., CX).
    * Column 3: ``[improper_flag]`` -- Set to 1 if the atom handles planar structures or improper dihedrals; set to 0 for sp3 non-planar centers.
    * Column 4: ``[group_id]`` -- Integer ID mapping the element group family (e.g., 0 for Carbon, 1 for Nitrogen, 2 for Oxygen, 4 for Hydrogen).
    * Column 5: ``[mass]`` -- Atomic floating-point mass (e.g., 12.01 or 14.01).
    * Column 6: ``[equivalent_flag]`` -- System index tracking fallback families. Match the header rules: use 1 for specialized types like cc/nc, 2 for types like cd/nd, or 0 for standard, typical environments.
    * Column 7: ``[atomic_num]`` -- Elemental integer atomic number (e.g., 6 for Carbon, 7 for Nitrogen).
* EQUA Line Fields:
    * Column 1: EQUA (Literal fallback string)
    * Column 2: ``[gaff_equivalent]`` -- The primary standard GAFF or AMBER atom type (e.g., c3, ca, h1) that matches your custom atom type's baseline structural parameters.
* CORR Line Fields (Optional but recommended for scoring):
    * Maps similarity index penalty scores when performing cross-field analogies. If you do not provide explicit CORR matrix lines for your new type, parmchk2 will fall back entirely onto the properties defined for your EQUA target.

The penalty matrix weights on the CORR lines represent how much a force field parameter (like a bond length or an angle) changes when you swap one atom type for another.
They are chemically derived mismatch scores that parmchk2 uses to calculate an empirical penalty when it has to guess missing parameters by substitution.


The 9 floating-point values specified on each ``CORR`` line in ``PARMCHK.DAT`` 
represent chemically derived mismatch scores. The ``parmchk2`` utility uses these 
weights to calculate empirical penalties when approximating missing force field 
parameters via substitution.

.. list-table:: CORR Parameter Weights Specification
   :widths: 10 30 60
   :header-rows: 1

   * - Position
     - Target Parameter
     - Property Measured
   * - 1
     - Bond Length Force Constant (:math:`K_r`)
     - Stiffness penalty encountered when swapping atoms in a bond.
   * - 2
     - Equilibrium Bond Length (:math:`r_0`)
     - Physical length structural penalty for atom substitution in a bond.
   * - 3
     - Bond Angle Force Constant (:math:`K_\theta`)
     - Stiffness penalty when swapping the *outer* atoms of a 3-atom angle.
   * - 4
     - Bond Angle Force Constant (:math:`K_\theta`)
     - Stiffness penalty when swapping the *center* atom of a 3-atom angle.
   * - 5
     - Equilibrium Bond Angle (:math:`\theta_0`)
     - Angle degree variance penalty when swapping the *outer* atoms.
   * - 6
     - Equilibrium Bond Angle (:math:`\theta_0`)
     - Angle degree variance penalty when swapping the *center* atom.
   * - 7
     - Dihedral Torsion Barrier (:math:`V_n`)
     - Rotational energy barrier penalty for a 4-atom torsion twist.
   * - 8
     - Improper Dihedral Torsion
     - Structural energy penalty for out-of-plane planar deviations.
   * - 9
     - Cumulative Fit Index
     - The total aggregated baseline difference between the atom types.

.. note::
   A weight value of ``-1.0`` signifies that insufficient statistical data was available 
   in the source training sets to calculate an exact substitution variance. In these instances, 
   ``parmchk2`` safely defaults to a built-in backup constant.
   

Example:

 .. code-block::

    PARM    C       1               0       12.01   0       6
    EQUA    c
    CORR    c2   4.9  21.9   4.3   2.9   2.9   2.3  96.7  -1.0  36.1
    CORR    ca   3.1  13.8   3.8   2.9   3.0   2.0 109.9  -1.0  36.6
    CORR    cp   6.2  27.6   0.2   6.5   3.1   2.5 113.5  -1.0  40.1
    CORR    cq   6.2  27.6   0.2   6.5   3.1   2.5 113.5  -1.0  40.1
    CORR    cc   3.2  14.2   5.5   2.5   3.7   2.3  48.2  -1.0  26.6
    CORR    cd   3.2  14.2   5.5   2.5   3.7   2.3  48.2  -1.0  26.6
    CORR    ce   3.6  16.1   4.6   1.8   3.1   2.4 126.9  -1.0  39.9
    CORR    cf   3.6  16.1   4.6   1.8   3.1   2.4 126.9  -1.0  39.9
    #
    PARM    CA      1               0       12.01   0       6
    EQUA    ca 
    CORR    CB   1.2   4.9   4.7   0.1   7.0   0.1  19.9  -1.0  19.8
    CORR    CC   0.3   0.6  -1.0  -1.0   0.8   0.1   3.4  -1.0  17.8
    CORR    CD   2.0   8.3  -1.0  -1.0   1.4   0.1  44.6  -1.0  26.3
    CORR    CK  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  37.8
    CORR    CM   1.8   8.0   1.9   0.1   2.1   0.1  58.0  -1.0  26.1
    CORR    CN   0.4   1.7   2.2   0.1   6.4   0.1   5.7  -1.0  16.5
    

The values in that matrix are not simple chemical indicators (like electronegativity or size); they are the statistical root-mean-square errors (RMSE) or percent absolute errors resulting from substituting one atom type for another across a massive validation training set of thousands of organic molecules. If you want to mathematically calculate your own CORR matrix weights for a custom atom type, you have to execute a formal, script-driven force field optimization pipeline.

.. _calculating_corr_matrix:

==================================================
Calculating Custom CORR Penalty Matrix Weights
==================================================

The numerical values in the ``CORR`` penalty matrix represent statistical 
variance metrics—specifically root-mean-square errors (RMSE) or percent absolute 
errors—derived from substituting one atom type for another. 

To mathematically calculate and calibrate these weights for a novel custom atom 
type, you must execute a formal, data-driven force field parametrization 
pipeline.

.. contents:: Section Contents
   :local:
   :depth: 2

---------------------------------------------
Step 1: Assemble a Diverse Training Database
---------------------------------------------

A penalty matrix value cannot be computed from a single molecule. You must 
construct a diverse database consisting of multiple distinct chemical structures 
(typically 50 to 100 configurations) that incorporate your custom atom type 
across a variety of neighboring chemical environments.

-------------------------------------------
Step 2: Generate Quantum Mechanical Data
-------------------------------------------

For every molecule in your training database, you must generate a highly precise 
Quantum Mechanical (QM) reference profile:

1. Perform a high-level geometry optimization using a quantum engine (e.g., ORCA, Gaussian, or Q-Chem).
2. Execute a rigid or relaxed Potential Energy Surface (PES) scan across the targeted structural parameters (e.g., stretching a bond, bending an angle, or systematically rotating a dihedral angle).
3. Extract the exact QM energies (:math:`E_{\text{QM}}`) and force vectors (:math:`F_{\text{QM}}`) corresponding to each step along the coordinate scan.

--------------------------------------------------
Step 3: Evaluate Molecular Mechanics Substitutions
--------------------------------------------------

Replicate the exact configurations generated during the QM scans within a 
classical Molecular Mechanics (MM) framework:

1. Create a baseline topology where your custom atom type is explicitly replaced by the standard fallback atom type (e.g., substituting your custom type with standard ``c3``).
2. Evaluate the classical MM energy (:math:`E_{\text{MM}}`) and analytical forces (:math:`F_{\text{MM}}`) for each saved geometry step using a simulation engine (such as AMBER's ``sander`` module). 

.. warning::
   The structural coordinates must remain completely rigid and identical to the QM configurations during this single-point evaluation phase.

-------------------------------------------
Step 4: Compute Objective Deviation Errors
-------------------------------------------

The values for each column of the ``CORR`` matrix line are populated by 
quantifying the explicit deviation between the substitute MM model and the true 
QM baseline profiles.

To calculate the structural and stiffness error limits (such as the Bond Length 
Force Constant and Equilibrium Length errors at positions 1 and 2), apply an 
:math:`L_2` regularization algorithm to minimize the error objective function:

.. math::

   \chi^2 = \frac{1}{N} \sum_{i=1}^{N} \left( E_{\text{MM}}(\theta) - E_{\text{QM}} \right)^2

Where :math:`\theta` encapsulates the force constant parameters (e.g., :math:`K_r`) and structural equilibria (e.g., :math:`r_0`).

* **Equilibrium Deviations (Positions 2, 5, 6):** Calculated as the absolute physical difference or percent error between the QM-optimized minimum structural value and the substitute baseline parameter.
* **Force Constant Curvature Error (Positions 1, 3, 4, 7, 8):** Quantifies the explicit variance in the parabolic curvature of the energy profiles when overlaying the QM scan against the substitute classical MM evaluation.

--------------------------------------------
Step 5: Automate Calibrations via paramfit
--------------------------------------------

Rather than computing these mathematical error terms manually from scratch, it 
is highly recommended to feed your structural files and QM energy logs 
directly into AMBER's native parameter-fitting module, ``paramfit``.

.. code-block:: bash

   # Example invocation structure for parameter fitting
   paramfit -i paramfit.in -p system.prmtop -c system.inpcrd > paramfit.out

By providing your targeted QM energy profiles alongside the alternative system 
topologies, ``paramfit`` automatically runs the sum-of-squares evaluations. 
The resulting statistical errors can then be normalized and pasted directly 
into your custom ``CORR`` rows within your local ``PARMCHK.DAT`` file.

    
    
Tutorial
--------


1. User options
~~~~~~~~~~~~~~~

 .. code-block::

    usage: ffpopt-parmchk2.py [-h] -i INP -o OUT [-s {1,gaff,2,gaff2,3,parm99,4,parm10,5,lipid14}] [--ffmiss {1,gaff,2,gaff2}]
                          [-atc SCORE_FILE] [-frc PARMCHK_FRC] [-afrc PARMCHK_AFRC] [-a] [-v] [--clean] [--only-missing]

    ffpopt parameterization wrapper mapping traditional parmchk2 CLI flags to RunCreateAmberFrcmod.

    options:
       -h, --help            show this help message and exit
       -i INP, --inp INP, --input-file INP
                             Path to the molecular structure input target file or coordinate string context.
       -o OUT, --out OUT, --output-file OUT
                             Path location where the generated output parameter modification (.frcmod) file will be saved.
       -s {1,gaff,2,gaff2,3,parm99,4,parm10,5,lipid14}, --fftype {1,gaff,2,gaff2,3,parm99,4,parm10,5,lipid14}
                             Baseline Amber force field parameters selection template matched to your atom types. Complete valid
                             choices: 1, gaff, 2, gaff2, 3, parm99, 4, parm10, 5, lipid14. Default: gaff.
       --ffmiss {1,gaff,2,gaff2}
                             Baseline Amber force field used to replace missing ATTN parameters. Complete valid choices: 1, gaff, 2,
                             gaff2. Default: gaff.
       -atc SCORE_FILE, --score-file SCORE_FILE
                             Additional atom type score file, optional. Type 'parmchk2 -l' for details
       -frc PARMCHK_FRC, --ffstd PARMCHK_FRC
                             Optional string containing user-supplied custom frcmod/parameter
			     template paths. Multiple entries can be chained via a '+' delimiter
			     (e.g., 'ff14SB+ol15'). The names are case-sensitive and parmchk2
			     will silently drop invalid names. The valid names are: ff14SB,
			     ff99SB, ff03, bsc1, ol15, yil. The ol3 force field uses the same
			     frcmod as ol15. The base force field for each of these is parm10,
			     except ff99SB and f03, which use parm99; however, parm99 is
       -afrc PARMCHK_AFRC, --ffuser PARMCHK_AFRC
                             Optional input frcmod file used to supplement the parameters within -s and -frc.
       -a, --all-parameters  If present, prints all force field parameters for the system to the output file instead of only missing
                             records.
       -v, -vv, --verbose    Turn on explicit verbose console tracing execution outputs.
       --clean               Purge individual intermediate runtime temporary sandbox environments and report files when finished.
       --only-missing        Do not insert MASS and NONBON entries for those atoms within the standard force fields. This is the
                             default behavior of parmchk2; however, parmed will likely fail to manipulate and write the frcmod without
                             these entries.



			     
2. Execution Examples
~~~~~~~~~~~~~~~~~~~~~

The following examples will use parmchk2 to create a frcmod file containing
the force field parameters that are not within the base and standard force fields.
It will then use the gaff force field to fill-in any parameters that parmchk2
marked with "ATTN".
To make the output frcmod file self-contained in a manner that parmed can
read, manipulate, and write the frcmod file, the frcmod file must contain
entries in the MASS and NONBON sections for all types appearing anywhere within
the frcmod file.
The default behavior of ffpopt-parmchk2.py is to include the entries in MASS
and NONBON.

a. Basic parameter discovery using short options:
   
 .. code-block::
     
   [user@code] ffpopt-parmchk2.py -i ligand.mol2 -o ligand.frcmod --fftype gaff2

b. Complex parameter assembly using long options and custom force fields string:

 .. code-block::
      
   [user@code] ffpopt-parmchk2.py --inp non_canonical_rna.mol2 \
        --out system.frcmod \
        --fftype "parm10" \
        --ffstd "ff14SB+ol15+ol3" \
	--ffmiss "gaff" \
        --verbose \
        --clean
	

.. note::
   If you want to reproduce the default behavior of parmchk2, which
   only includes MASS and NONBON entries if they were missing in the
   base and standard force fields, then use the ``--only-missing`` option.

   
   
