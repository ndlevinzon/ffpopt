.. _tleap-tutorial:


Preparing lib, parm7, and rst7 files with ffpopt-tleap.py
=========================================================

This tutorial will show how to prepare a parm7, rst7, and lib
file for amber using ffpopt-tleap.py.
This script is also capable of preparing a tar.gz file containing
the tleap script and all necessary inputs to prepare the files
from scratch.


Learning Objectives
-------------------


Tutorial
--------


1. User Options
~~~~~~~~~~~~~~~

 .. code-block::

    usage: ffpopt-tleap.py [-h] -i INP [-s {1,gaff,2,gaff2,3,parm99,4,parm10,5,lipid14}]
           [-frc PARMCHK_FRC] [--frcmod FRCMOD_LIST] [--tar] [--verbose] [--clean]

    ffpopt system orchestration wrapper mapping CLI configurations to RunBuildAmberSystem.

    options:
      -h, --help            show this help message and exit
      -i INP, --inp INP, --input-file INP
                            Path to the molecular structure input target (e.g., a file path
			    string or coordinate state context).
      -s {1,gaff,2,gaff2,3,parm99,4,parm10,5,lipid14},
      --fftype {1,gaff,2,gaff2,3,parm99,4,parm10,5,lipid14}
                            Baseline Amber force field parameters selection template matched
			    to your atom types. Complete valid choices: 1, gaff, 2, gaff2,
			    3, parm99, 4, parm10, 5, lipid14. Default: gaff.
      -frc PARMCHK_FRC, --ffstd PARMCHK_FRC
                            Optional string containing user-supplied custom frcmod/parameter
			    template paths. Multiple entries can be chained via a '+' delimiter
			    (e.g., 'ff14SB+ol15'). The names are case-sensitive and parmchk2
			    will silently drop invalid names. The valid names read by parmchk2
			    are: ff14SB, ff99SB, ff03, bsc1, ol15, yil. The ol3 force field
			    uses the same frcmod as ol15. The base force field for each of
			    these is parm10, except ff99SB and f03, which use parm99; however,
			    parm99 is deprecated.  The valid names used by tleap are:
			    ff99SB (deprecated), ff03 (deprecated), ff14SB, fb15, ff19SB,
			    bsc1, ol15, ol21, ol24, ol3, roc, shaw, LJbb, yil.
      -afrc, --frcmod FRCMOD_LIST
                            Explicit file path to a pre-generated .frcmod modification file.
                            This option can be used 0-or-more times to inject multiple files
			    sequentially into LEaP.
      --tar                 If present, bundles all successfully generated system physics
                            parameters, topology files, and scripts into a consolidated
			    archive.
      -v, --verbose         Turn on explicit verbose console tracing execution outputs.
      --clean               Purge individual intermediate runtime temporary sandbox
                            environments when finished.


The ``--fftype`` and ``--ffstd`` options refer to standard Amber force fields. There are
the same options used in ffpopt-parmchk2.py. In the present context, these are options
are used to control what amber force field files to load.

The ``--frc`` (or ``--ffstd``) option plays a duel role. Some of the settings are passed
directly to parmchk2 (but this only happens if ``--frcmod`` (or ``-afrc``) is not used.
The ``--frc`` option is also use to source leaprc inputs when running tleap.

.. list-table:: -\-ffstd options
   :widths: 20 40 40
   :header-rows: 1

   * - -\-ffstd name
     - Effect in tleap
     - Effect in parmchk2
   * - ff99SB
     - source leaprc.protein.ff99SB
     - imports dat/leap/parm/frcmod.ff99SB
   * - ff03
     - source leaprc.protein.ff03.r1
     - imports dat/leap/parm/frcmod.ff03
   * - ff14SB
     - source leaprc.protein.ff14SB
     - imports dat/leap/parm/frcmod.ff14SB
   * - fb15
     - source leaprc.protein.fb15
     - 
   * - ff19SB
     - source leaprc.protein.ff19SB
     - 
   * - bsc1
     - source leaprc.DNA.bsc1
     - imports dat/leap/parm/frcmod.parmbsc1
   * - ol15
     - source leaprc.DNA.OL15
     - imports dat/leap/parm/frcmod.DNA.OL15
   * - ol21
     - source leaprc.DNA.OL21
     - 
   * - ol24
     - source leaprc.DNA.OL24
     - 
   * - ol3
     - source leaprc.RNA.OL3
     - (use ol15+ol3 instead of just ol3)
   * - roc
     - source leaprc.RNA.ROC
     - 
   * - shaw
     - source leaprc.RNA.Shaw
     - 
   * - LJbb
     - source leaprc.RNA.LJbb
     - 
   * - yil
     - source leaprc.RNA.YIL
     - imports dat/leap/parm/frcmod.parmCHI_YIL



2. Basic automated pipeline (generates .frcmod, .lib, .parm7, .rst7 via GAFF2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mol.json file contains conformers of a residue named MOL. The first conformer
is extracted as a mol2 file and run through tleap. In this usage, a separate
MOL.frcmod file is not specified, so a MOL.frcmod is created via internal library
calls to ffpopt-parmchk2.py with the provided force field options.

 .. code-block::

    [user@comp] ffpopt-tleap.py --inp mol.json --fftype gaff2
    [user@comp] ls MOL*
    MOL.frcmod  MOL.lib  MOL.mol2  MOL.parm7  MOL.rst7  MOL.sh
    [user@comp] ffpopt-tleap.py --inp mol.json --fftype gaff2 --tar
    [user@comp] ls MOL*
    MOL.frcmod  MOL.lib  MOL.mol2  MOL.parm7  MOL.rst7  MOL.sh  MOL.tar.gz
    [user@comp] mkdir tmp; cd tmp; tar -xzf ../MOL.tar.gz; ls
    MOL.frcmod  MOL.lib  MOL.mol2  MOL.parm7  MOL.rst7  MOL.sh  tleap_commands.in


3. Loading multiple pre-generated .frcmod files explicitly (skips parmchk2 step)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 .. code-block::

    [user@comp] ffpopt-tleap.py -i mol.mol2 \
             --fftype parm10 \
	     --ffstd ff14SB+ol3 \
	     --frcmod modification_A.frcmod \
	     --frcmod modification_B.frcmod \
	     --tar --verbose --clean


4. Update the json to use the parm7
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 .. code-block:
    
    [user@comp] ffpopt-PrepareInput.py --update --parm MOL.parm7 --crd mol.json --out mol.json


Now you can use mol.json with ``ffpopt-Optimize.py --model=sander --inp mol.json --out mol.sander.json``


5. Example workflow for building a molecule with GAFF2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 .. code-block::

    [user@comp] ffpopt-antechamber.py -i ${RES}.pdb:${CHARGE} -o ${RES}.mol2 -rn ${RES} \
         -c abcg2 -at gaff2 -v --clean
    [user@comp] ffpopt-tleap.py -i ${RES}.mol2 -s gaff2 --tar -v --clean 

    
6. Example workflow for building a molecule with ff14SB+OL15+OL3
   
 .. code-block::
    
    [user@comp] ffpopt-antechamber.py -i ${RES}.pdb:${CHARGE} -o ${RES}.mol2 -rn ${RES} \
         -c am1bcc -at amber -v --clean
    [user@comp] ffpopt-tleap.py -i ${RES}.mol2 -s parm10 -frc "ff14SB+ol15" --ffmiss gaff2 --tar -v --clean 


The above example is equivalent to the following

 .. code-block::
    
    [user@comp] ffpopt-antechamber.py -i ${RES}.pdb:${CHARGE} -o ${RES}.mol2 -rn ${RES} \
         -c am1bcc -at amber -v --clean
    [user@comp] ffpopt-parmchk2.py -i ${RES}.mol2 -o ${RES}.frcmod -s parm10 \
         -frc "ff14SB+ol15+ol3" --ffmiss gaff2 -v --clean
    [user@comp] ffpopt-tleap.py -i ${RES}.mol2 -s parm10 -frc "ff14SB+ol15+ol3" \
         -afrc ${RES}.frcmod --tar -v --clean 

