#!/usr/bin/env python3
"""Command-line utility wrapper for the ffpopt parmchk2 workflow.

This script provides a unified command-line interface (CLI) that translates traditional 
AmberTools 'parmchk2' syntax options into the ffpopt framework. It wraps the underlying
`RunCreateAmberFrcmod` library function to scan incoming structures for missing 
force field valence parameters and generate corresponding modification (.frcmod) files.
"""

import sys

def main():
    """Parse compatible parmchk2 command-line flags and execute RunCreateAmberFrcmod."""
    import argparse
    import traceback

    from ffpopt.CreateAmberLigand import RunCreateAmberFrcmod

    # Define a custom RawDescriptionHelpFormatter to preserve multi-line text blocks in help menus
    parser = argparse.ArgumentParser(
        description="Ffpopt parameterization wrapper mapping traditional parmchk2 CLI flags to RunCreateAmberFrcmod.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Execution Examples:
------------------
1. Basic parameter discovery using short options:
   $ ./ffpopt-parmchk2.py -i ligand.mol2 -o ligand.frcmod -at gaff2

2. Complex parameter assembly using long options and custom force fields string:
   $ ./ffpopt-parmchk2.py --inp non_canonical_rna.mol2 \\
                          --out system.frcmod \\
                          --fftype "parm10" \\
                          --ffstd "ff14SB+ol15" \\
                          --verbose \\
                          --clean

Output Generation Notes:
------------------------
This utility reads a molecular structure or file (such as a Mol2 or SMILES string),
extracts structural connectivity metrics, and executes parmchk2. It outputs a
Force Field Modification (.frcmod) file containing missing parameter definitions 
reconstructed with reasonable parameter suggestions.
"""
    )

    # Map input and output parameter flags, accepting short and dash-separated long variants
    parser.add_argument("-i", "--inp", "--input-file", required=True, dest="inp", 
                        help="Path to the molecular structure input target file or coordinate string context.")
    parser.add_argument("-o", "--out", "--output-file", required=True, dest="out", 
                        help="Path location where the generated output parameter modification (.frcmod) file will be saved.")

    # Commented out parameters not currently referenced by the underlying library target signatures
    # parser.add_argument("-f", "--input-format", default="mol2", help="Input file format (prepi, prepc, ac, mol2, frcmod, leaplog).")
    # parser.add_argument("-p", "--parmfile", help="Raw baseline parameter dat file path reference, overrides -at.")

    # Map optional parameter flags with clean matching long-option forms and complete value constraints
    parser.add_argument("-s", "--fftype", default="gaff", dest="parmchk_ff",
                        choices=["1", "gaff", "2", "gaff2", "3", "parm99", "4", "parm10", "5", "lipid14"],
                        help="Baseline Amber force field parameters selection template matched to your atom types. Complete valid choices: 1, gaff, 2, gaff2, 3, parm99, 4, parm10, 5, lipid14. Default: gaff.")

    parser.add_argument("--ffmiss", default="gaff", dest="parmchk_ffmissing",
                        choices=["1", "gaff", "2", "gaff2"],
                        help="Baseline Amber force field used to replace missing ATTN parameters. Complete valid choices: 1, gaff, 2, gaff2. Default: gaff.")


    parser.add_argument("-atc", "--score-file",
                        required=False,
                        help="Additional atom type score file, optional. Type 'parmchk2 -l' for details")

    
    parser.add_argument("-frc", "--ffstd", dest="parmchk_frc",
                        help="Optional string containing user-supplied custom frcmod/parameter template paths. Multiple entries can be chained via a '+' delimiter (e.g., 'ff14SB+ol15'). The names are case-sensitive and parmchk2 will silently drop invalid names. The valid names are: ff14SB, ff99SB, ff03, bsc1, ol15, yil.  The ol3 force field uses the same frcmod as ol15. The base force field for each of these is parm10, except ff99SB and f03, which use parm99; however, parm99 is deprecated.")

    parser.add_argument("-afrc", "--ffuser", dest="parmchk_afrc",
                        help="Optional input frcmod file used to supplement the parameters within -s and -frc.")

    
    
    parser.add_argument("-a", "--all-parameters",
                        action='store_true',
                        dest="parmchk_all",
                        help="If present, prints all force field parameters for the system to the output file instead of only missing records.")

    parser.add_argument("-v","-vv","--verbose",
                        action='store_true',
                        help="Turn on explicit verbose console tracing execution outputs.")
    
    parser.add_argument("--clean",
                        action='store_true',
                        help="Purge individual intermediate runtime temporary sandbox environments and report files when finished.")

    parser.add_argument("--only-missing",
                        action='store_true',
                        help="Do not insert MASS and NONBON entries for those atoms within the standard force fields. This is the default behavior of parmchk2; however, parmed will likely fail to manipulate and write the frcmod without these entries.")


    
    # Parse args array 
    args = parser.parse_args()

    verbose_bool = args.verbose
    clean_bool = args.clean

    if verbose_bool:
        print(f"[ffpopt-parmchk2] Launching parameters validation script loop targeting: {args.inp}")

    try:
        # Call the underlying core library function
        res = RunCreateAmberFrcmod(
            inp=args.inp,
            out=args.out,
            clean=clean_bool,
            verbose=verbose_bool,
            parmchk_ff=args.parmchk_ff,
            parmchk_frc=args.parmchk_frc,
            parmchk_keep=not clean_bool,  # Evaluated seamlessly directly from clean state status
            parmchk_all=args.parmchk_all,
            parmchk_ffmissing=args.parmchk_ffmissing,
            parmchk_atc=args.score_file,
            parmchk_afrc=args.parmchk_afrc,
            only_missing=args.only_missing
        )

        if verbose_bool:
            print(f"[ffpopt-parmchk2] Processing completed successfully. Output parameter saved to: {args.out}")

    except Exception as e:
        # Guarantee trace tracking backplanes output directly regardless of hidden status parameters flags
        print(f"[ffpopt-parmchk2] Critical Error: Underlying parameters discovery loop failed. Reason: {str(e)}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
