#!/usr/bin/env python3
"""Command-line utility wrapper for the ffpopt antechamber workflow.

This script provides a unified command-line interface (CLI) that translates traditional 
AmberTools 'antechamber' syntax options into the ffpopt framework. It handles single 
or multiple conformation sets, executes small molecule charge definitions, and performs 
conformer-averaged calculations by wrapping the underlying `RunCreateAmberMol2` library function.
"""

import sys

def main():
    """Parse compatible antechamber command-line flags and execute RunCreateAmberMol2."""
    import argparse
    import traceback

    from ffpopt.CreateAmberLigand import RunCreateAmberMol2

    # Define a custom RawDescriptionHelpFormatter to preserve multi-line text blocks in help menus
    parser = argparse.ArgumentParser(
        description="Ffpopt parameterization wrapper mapping traditional antechamber CLI flags to RunCreateAmberMol2.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Execution Examples:
------------------
1. Basic small molecule parameterization using short options:
   $ ./ffpopt-antechamber.py -i ligand.xyz -o ligand.mol2 -rn LIG -at gaff2 -c abcg2

2. Comprehensive conformer-averaged parameterization using long options:
   $ ./ffpopt-antechamber.py --inp molecule.mol2 \\
                             --out processed.mol2 \\
                             --resname MOL \\
                             --atom-type gaff2 \\
                             --charge-method am1bcc \\
                             --nconf 100 \\
                             --nkeep 20 \\
                             --verbose --clean \\

Output Generation Notes:
------------------------
This utility creates two sets of file paths:
  - The main file path targeting your designated output destination ([output_file].mol2) 
    retaining conformer-averaged atomic charges.
  - A sister file path context generated on disk appended with '.confdepqs.mol2' 
    retaining unique, conformer-dependent raw charge matrices.
"""
    )

    # Map input and output parameter flags, accepting short and dash-separated long variants
    parser.add_argument("-i", "--inp", "--input-file", required=True, dest="inp", 
                        help="Path to the molecular structure input target (e.g., .xyz, .pdb, .mol2, .json, or inchi/smiles). If the input is an xyz or pdb, then the name can be suffixed with a charge. Example: '--inp=foo.xyz:-1' means the net charge of the molecule is -1. Example: '--inp=foo.pdb:1' means the net charge of the molecule is +1. The default is to assume the molecule is neutral.")
    parser.add_argument("-o", "--out", "--output-file", required=True, dest="out", 
                        help="Output .mol2 or .json")

    # Commented out parameters not currently referenced by the underlying library target signatures
    # parser.add_argument("-fi", "--input-format", required=True, help="Input file format descriptor.")
    # parser.add_argument("-fo", "--output-format", required=True, help="Output file format descriptor.")
    # parser.add_argument("-nc", "--net-charge", default="0", help="Net molecular charge (int).")
    # parser.add_argument("-m", "--multiplicity", default="1", help="Multiplicity (2S+1).")

    # Map optional parameter flags with clean matching long-option forms and complete value constraints
    parser.add_argument("-c", "--charge-method", default="am1bcc", dest="charge_method",
                        choices=["bcc", "am1bcc", "abcg2"],
                        help="Antechamber charge evaluation routing scheme choice. Complete valid choices: bcc, am1bcc, abcg2. Default: am1bcc.")
    
    parser.add_argument("-at", "--atom-type", default="gaff", dest="atom_type",
                        choices=["gaff", "gaff2", "amber", "bcc", "abcg2", "sybyl"],
                        help="Force field atom typing family classification matrix rules. Complete valid choices: gaff, gaff2, amber, bcc, abcg2, sybyl. Default: gaff.")
    
    parser.add_argument("-rn", "--resname", default="MOL", dest="residue_name", 
                        help="Residue designator code naming convention token to label the system outputs. Default: MOL.")
    
    parser.add_argument("-v","--verbose",
                        action='store_true',
                        help="Turn on verbose printing.")
    
    parser.add_argument("--clean",
                        action='store_true',
                        help="Purge intermediate runtime file wrappers and scratch environments when finished.")

    # Custom ffpopt extensions supporting the embedded conformer search steps
    parser.add_argument("--nconf", type=int, default=0, 
                        help="Number of initial candidate configuration geometries to generate inside the conformer pool search loop. Set to 0 to skip entirely. Default: 0")
    parser.add_argument("--nkeep", type=int, default=100, 
                        help="Maximum structural tracking storage threshold count used to clamp the kept lowest-energy ensemble pool members. Default: 100")
    
    parser.add_argument("--no-opt",
                        action='store_true',
                        help="If present, don't perform a am1 optimization when calculating charges")

    parser.add_argument("--confdep",
                        action='store_true',
                        help="If present, save conformer-dependent charges to {base}.confdepqs.{ext} where base and ext are the base name and extension of --out")

    # Parse args array 
    args = parser.parse_args()

    # Translate your parsed CLI options variables directly into your library's unified keyword dictionary entries
    # -s / --verbose -> verbose (True/False based on level choice strings)
    # -pf / --clean -> clean (y maps to True, n maps to False)
    #verbose_bool = True if args.s in ["1", "2"] else False
    #clean_bool = True if args.pf.lower() == "y" else False
    verbose_bool = args.verbose
    clean_bool = args.clean

    if verbose_bool:
        print(f"[ffpopt-antechamber] Launching utility wrapper tool sequence pipeline path targeting: {args.inp}")

    try:
        # Call the underlying core function that you supplied
        avglos, los = RunCreateAmberMol2(
            inp=args.inp,
            out=args.out,
            clean=clean_bool,
            verbose=verbose_bool,
            confsearch_nconf=args.nconf,
            confsearch_nkeep=args.nkeep,
            antechamber_resname=args.residue_name,
            antechamber_atomtype=args.atom_type,
            antechamber_chargetype=args.charge_method,
            antechamber_optimize=not args.no_opt,
            confdep=args.confdep
        )

        if verbose_bool:
            print(f"[ffpopt-antechamber] Processing completed successfully. Output files saved via base name track context: {args.out}")

    except Exception as e:
        # Guarantee trace tracking backplanes output directly regardless of hidden status parameters flags
        print(f"[ffpopt-antechamber] Critical Error: Underlying execution wrapper failed. Reason: {str(e)}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
