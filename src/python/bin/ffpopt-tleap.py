#!/usr/bin/env python3
"""Command-line utility wrapper for the ffpopt system orchestration workflow.

This script provides a unified command-line interface (CLI) that translates traditional 
AmberTools parameters into the ffpopt automation framework. It wraps the high-level
`RunBuildAmberSystem` function to coordinate parmchk2 parameter validation and tleap 
topology creation, and supports passing pre-computed .frcmod files 0-or-more times.
"""

import sys

def main():
    """Parse compatible parameters and execute RunBuildAmberSystem orchestrator."""
    import argparse
    import traceback

    from ffpopt.CreateAmberLigand import RunBuildAmberSystem

    # Define a custom RawDescriptionHelpFormatter to preserve multi-line text blocks in help menus
    parser = argparse.ArgumentParser(
        description="Ffpopt system orchestration wrapper mapping CLI configurations to RunBuildAmberSystem.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Execution Examples:
------------------
1. Basic automated pipeline (generates .frcmod, .lib, .parm7, .rst7 via GAFF2):
   $ ./ffpopt-tleap.py --inp ligand.mol2 --fftype gaff2

2. Loading multiple pre-generated .frcmod files explicitly (skips parmchk2 step):
   $ ./ffpopt-tleap.py -i core.mol2 \\
                       --fftype parm10 \\
                       --ffstd ff14SB+ol3 \\
                       --frcmod modification_A.frcmod \\
                       --frcmod modification_B.frcmod \\
                       --tar \\
                       --verbose \\
                       --clean

Output Generation Notes:
------------------------
This utility creates an uncompressed set of simulation-ready files (.mol2, .frcmod, 
.lib, .parm7, .rst7, .sh, and tleap_commands.in) directly inside the current directory. 
If the '--tar' flag is active, a consolidated and flattened '[resname].tar.gz' 
archive is generated alongside them.
"""
    )

    # Map the primary mandatory structural input parameter flag, accepting both short and long variants
    parser.add_argument("-i", "--inp", "--input-file", required=True, dest="inp", 
                        help="Path to the molecular structure input target (e.g., a file path string or coordinate state context).")

    # Map optional parameter flags with clean matching long-option forms and complete value constraints
    parser.add_argument("-s", "--fftype", default="gaff", dest="parmchk_ff",
                        choices=["1", "gaff", "2", "gaff2", "3", "parm99", "4", "parm10", "5", "lipid14"],
                        help="Baseline Amber force field parameters selection template matched to your atom types. Complete valid choices: 1, gaff, 2, gaff2, 3, parm99, 4, parm10, 5, lipid14. Default: gaff.")
    
    parser.add_argument("-frc", "--ffstd", dest="parmchk_frc",
                        help="Optional string containing user-supplied custom frcmod/parameter template paths. Multiple entries can be chained via a '+' delimiter (e.g., 'ff14SB+ol15'). The names are case-sensitive and parmchk2 will silently drop invalid names. The valid names are: ff99SB (deprecated), ff03 (deprecated), ff14SB, fb15, ff19SB, bsc1, ol15, ol21, ol24, ol3, roc, shaw, LJbb, yil.")

    # Implement the 0-or-more array append structure using action='append'
    parser.add_argument("-afrc", "--frcmod", action='append', dest="frcmod_list",
                        help="Explicit file path to a pre-generated .frcmod modification file. This option can be used 0-or-more times to inject multiple files sequentially into LEaP.")

    parser.add_argument("--ffmiss", default="gaff", dest="parmchk_ffmissing",
                        choices=["1", "gaff", "2", "gaff2"],
                        help="Baseline Amber force field used to replace missing ATTN parameters. This is only used if --frcmod/-afrc is not used. Complete valid choices: 1, gaff, 2, gaff2. Default: gaff.")

    
    parser.add_argument("--tar",
                        action='store_true',
                        help="If present, bundles all successfully generated system physics parameters, topology files, and scripts into a consolidated archive.")

    parser.add_argument("-v","--verbose",
                        action='store_true',
                        help="Turn on explicit verbose console tracing execution outputs.")
    
    parser.add_argument("--clean",
                        action='store_true',
                        help="Purge individual intermediate runtime temporary sandbox environments when finished.")

    # Parse args array 
    args = parser.parse_args()

    verbose_bool = args.verbose
    clean_bool = args.clean
    tar_bool = args.tar

    if verbose_bool:
        print(f"[ffpopt-tleap] Launching system topology builder orchestrator sequence targeting: {args.inp}")

    try:
        # Pass variables cleanly straight into your high-level orchestrator function block
        topology_results = RunBuildAmberSystem(
            inp=args.inp,
            parmchk_ff=args.parmchk_ff,
            parmchk_frc=args.parmchk_frc,
            parmchk_ffmissing=args.parmchk_ffmissing,
            frcmod_files=args.frcmod_list,  # Passes None or the collected array list of file strings
            clean=clean_bool,
            verbose=verbose_bool,
            tar=tar_bool
        )

        if verbose_bool:
            print("[ffpopt-tleap] Pipeline completed successfully. Output files saved in the current directory.")
            if "tar" in topology_results:
                print(f"[ffpopt-tleap] Consolidated archive bundle published to: {topology_results['tar']}")

    except Exception as e:
        # Guarantee trace tracking backplanes output directly regardless of hidden status parameters flags
        print(f"[ffpopt-tleap] Critical Error: System builder orchestration loop encountered an error. Reason: {str(e)}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
