#!/usr/bin/env python3

import os
import subprocess
import tempfile
import shutil
import copy
from pathlib import Path
from types import SimpleNamespace

# Central authority matching parmchk2 help definitions to modern tleap macro scripts
_BASE_FF_MAP = {
    "1":        {"parmchk": "1",        "tleap": ["source leaprc.gaff"]},
    "gaff":     {"parmchk": "1",     "tleap": ["source leaprc.gaff"]},
    "2":        {"parmchk": "2",        "tleap": ["source leaprc.gaff2"]},
    "gaff2":    {"parmchk": "2",    "tleap": ["source leaprc.gaff2"]},
    "3":        {"parmchk": "3",        "tleap": ["loadAmberParams parm99.dat"]},
    "parm99":   {"parmchk": "3",   "tleap": ["loadAmberParams parm99.dat"]},
    "4":        {"parmchk": "4",        "tleap": ["loadAmberParams parm10.dat"]},
    "parm10":   {"parmchk": "4",   "tleap": ["loadAmberParams parm10.dat"]},
    "5":        {"parmchk": "5",        "tleap": ["source leaprc.lipid14"]},
    "lipid14":  {"parmchk": "5",  "tleap": ["source leaprc.lipid14"]},
}

_MODIFIER_FF_MAP = {
    # Proteins
    "ff14sb":   {"tleap": "source leaprc.protein.ff14SB"},
    "ff99sb":   {"tleap": "source leaprc.protein.ff99SB"},
    "ff03":     {"tleap": "source leaprc.protein.ff03.r1"},
    "fb15":     {"tleap": "source leaprc.protein.fb15"},
    "ff19sb":   {"tleap": "source leaprc.protein.ff19SB"},
    # DNA
    "bsc1":     {"tleap": "source leaprc.DNA.bsc1"},
    "ol15":     {"tleap": "source leaprc.DNA.OL15"},
    "ol21":     {"tleap": "source leaprc.DNA.OL21"},
    "ol24":     {"tleap": "source leaprc.DNA.OL24"},
    # RNA
    "ol3":      {"tleap": "source leaprc.RNA.OL3"},
    "roc":      {"tleap": "source leaprc.RNA.ROC"},
    "shaw":     {"tleap": "source leaprc.RNA.Shaw"},
    "ljbb":     {"tleap": "source leaprc.RNA.LJbb"},
    "yil":      {"tleap": "source leaprc.RNA.YIL"},
}

def resolve_tleap_commands(parmchk_ff: str, parmchk_frc: str = None) -> list:
    """Helper tool to convert unified parmchk2 strings into explicit tleap script lines."""
    commands = []
    
    # Process Base Force Field
    clean_base = str(parmchk_ff).strip().lower()
    if clean_base in _BASE_FF_MAP:
        commands.extend(_BASE_FF_MAP[clean_base]["tleap"])
    else:
        # Fallback if user passes direct leaprc paths or custom dat parameters
        if "leaprc" in clean_base:
            commands.append(f"source {parmchk_ff}")
        else:
            commands.append(f"loadAmberParams {parmchk_ff}")

    # Process Extra Modification Tokens (e.g. "ff14SB+bsc1")
    if parmchk_frc:
        tokens = [t.strip() for t in parmchk_frc.split("+") if t.strip()]
        for token in tokens:
            token_lower = token.lower()
            if token_lower in _MODIFIER_FF_MAP:
                commands.append(_MODIFIER_FF_MAP[token_lower]["tleap"])
            else:
                if "leaprc" in token_lower:
                    commands.append(f"source {token}")
                else:
                    commands.append(f"loadAmberParams {token}")
                    
    return commands




def RunCreateAmberMol2\
        (*,
         inp: str,
         out: str,
         clean: bool = False,
         verbose: bool = True,
         confsearch_nconf = 0,
         confsearch_nkeep = 100,
         antechamber_resname: str = None,
         antechamber_atomtype: str = "gaff",
         antechamber_chargetype: str = "am1bcc",
         antechamber_optimize: bool = False,
         confdep: bool = False,
         **standard_kwargs ):
    """Generate Amber Mol2 files with conformer-averaged and conformer-dependent charges.

    This function reads a molecular structure or file, optionally performs a 
    conformer search, runs Antechamber to assign atom types and charges, 
    and saves the resulting outputs. It yields two sets of structures: one 
    with ensemble-averaged atomic charges across all conformers, and one 
    retaining individual conformer charges.

    Parameters
    ----------
    inp : str, Struct, or ListOfStruct
        Input molecular structure. Can be a filepath to a structural file, 
        a single `Struct` object, or a `ListOfStruct` collection.
    out : str
        Base filepath where the output conformer-averaged Mol2 file is saved.
        Conformer-dependent charges are saved to a separate file appended 
        with '.confdepqs'.
    clean : bool, default False
        If True, the temporary working directory is deleted upon completion.
        If False, files are left intact for inspection or debugging.
    verbose : bool, default True
        If True, prints progress details and temporary file paths to stdout.
    confsearch_nconf : int, default 0
        Number of initial conformers to generate during the conformer search.
        Set to 0 to skip the conformer search entirely.
    confsearch_nkeep : int, default 100
        Maximum number of low-energy conformers to keep from the search.
    antechamber_resname : str, optional (default None)
        Residue name assigned to the molecules in the output Mol2 files.
    antechamber_atomtype : str, default "gaff"
        Antechamber atom typing scheme to apply (e.g., 'gaff', 'gaff2').
    antechamber_chargetype : str, default "am1bcc"
        Antechamber charge assignment method (e.g., 'am1bcc', 'gasteiger').
    antechamber_optimize : bool, default False
        If True, instructs Antechamber to optimize molecular geometry 
        prior to quantum chemical charge calculation.
    confdep : bool, default False
        If True, save conformer-dependent charges to file named {base}.confdepqs.{ext}
        where pref and ext are the basename and extension of filename provided
        by "out". For example, if out="mol.json", then confdep=True will write
        mol.json with the conformer-averaged charges and mol.confdepqs.json with
        the conformer-dependent charges. One normally sets confdep=False.
    **standard_kwargs : dict
        Additional keyword arguments evaluated against standard internal model 
        defaults. Unexpected parameters raise a TypeError.

    Returns
    -------
    avglos : ListOfStruct
        A copy of the input structures containing ensemble conformer-averaged 
        atomic charges.
    los : ListOfStruct
        A copy of the input structures retaining unique conformer-dependent 
        atomic charges.

    Raises
    ------
    TypeError
        If an unexpected keyword argument is supplied through `**standard_kwargs`.

    See Also
    --------
    Antechamber.run_antechamber : Underlying charge and typing execution engine.
    confsearch.ConfSearch.ConformerSearch : Geometrical conformer generation.

    Examples
    --------
    >>> from my_module import RunCreateAmberMol2
    >>> avg_mol, dep_mol = RunCreateAmberMol2(
    ...     inp="ethanol.xyz",
    ...     out="ethanol.mol2",
    ...     confsearch_nconf=50,
    ...     antechamber_chargetype="am1bcc"
    ... )
    """

    
    import argparse
    from types import SimpleNamespace

    _p = argparse.ArgumentParser(add_help=False)
    #AddModelOptions(_p)
    std_defaults = vars(_p.parse_args([]))
    unknown = set(standard_kwargs) - set(std_defaults)
    if unknown:
        raise TypeError(
            f"Unexpected keyword argument(s): {sorted(unknown)}"
        )
    std = {**std_defaults, **standard_kwargs}

    args = SimpleNamespace(
        inp=inp,
        out=out,
        clean=clean,
        verbose=verbose,
        confsearch_nconf = confsearch_nconf,
        confsearch_nkeep = confsearch_nkeep,
        antechamber_resname=antechamber_resname,
        antechamber_atomtype=antechamber_atomtype,
        antechamber_chargetype=antechamber_chargetype,
        antechamber_optimize=antechamber_optimize,
        confdep=confdep,
        **std,
    )

    #########################################################################

    from . Antechamber import run_antechamber, FixCharges
    from . confsearch.ConfSearch import ConformerSearch
    from . Struct import ListOfStruct, Struct
    import os, copy, tempfile
    import numpy as np
    from pathlib import Path

    if isinstance(args.inp,ListOfStruct):
        los = copy.deepcopy(args.inp)
    elif isinstance(inp,Struct):
        los = ListOfStruct( [args.inp] )
    else:
        los = ListOfStruct.from_file(args.inp)

    if args.confsearch_nconf > 0 and args.confsearch_nkeep:
        los = ConformerSearch\
            (los,None,args.confsearch_nconf,
             args.confsearch_nkeep,True,250,
             0.5,not args.verbose)
    
    tmpbase = "tmpdir"
    os.makedirs(tmpbase, exist_ok=True)
    
    # 1. Manually create the unique temporary directory path
    tmpdir_path = tempfile.mkdtemp(dir=tmpbase)

    if args.verbose:
        print(f"Created temporary directory: {tmpdir_path}")
    
    try:

        for s in los:
            
            crds  = s.get_positions()
            elems = s.data["elements"]
            netq  = s.GetCharge()
            resname = s.data["resnames"][0]
            if args.antechamber_resname is not None:
                resname = args.antechamber_resname

            res = run_antechamber\
                ( crds, elems, netq,
                  get_charges=True,
                  get_types=True,
                  get_names=True,
                  charge_method=args.antechamber_chargetype,
                  atom_type_family=args.antechamber_atomtype,
                  optimize_geometry=args.antechamber_optimize,
                  adjust_names="n",
                  base_dir=tmpdir_path,
                  resname=resname,
                  verbose=args.verbose)
        
            s.data["charges"] = copy.deepcopy(res["charges"])
            s.data["types"] = copy.deepcopy(res["types"])
            s.data["names"] = copy.deepcopy(res["names"])
            for key in s.data["resnames"]:
                s.data["resnames"][key] = resname

    finally:
        # 4. Clean up only if requested
        if clean:
            if args.verbose:
                print(f"\nCleaning up: Removing {tmpdir_path}")
            import shutil
            shutil.rmtree(tmpdir_path)
        elif args.verbose:
            print(f"\n'clean' is False. Leaving directory intact at: {tmpdir_path}")


                
    charges = []
    for s in los:
        charges.append( s.data["charges"] )
        #print(s.data["charges"])
    charges = FixCharges( np.mean( np.array(charges), axis=0 ).tolist() )
    avglos = copy.deepcopy(los)
    for s in avglos:
        s.data["charges"] = copy.deepcopy(charges)


    if args.out is not None:
        avglos.save(args.out)
        if args.confdep:
            opath = Path(args.out)
            oname = opath.with_suffix(f".confdepqs{opath.suffix}")
            los.save( str(oname) )

    return avglos,los



def RunCreateAmberFrcmod\
        (*,
         inp: str, # mol2 or smiles or listofstruct
         out: str, # frcmod
         clean: bool = False,
         verbose: bool = True,
         parmchk_ff: str = "gaff",
         parmchk_frc: str = None,
         parmchk_keep: bool = False,
         parmchk_all: bool = False,
         parmchk_atc: str = None,
         parmchk_ffmissing: str = "gaff",
         parmchk_afrc: str = None,
         only_missing: bool = False,
         **standard_kwargs ):
    """Generate an Amber frcmod file using the parmchk2 utility.

    This function reads a molecular structure or file (such as a Mol2 or 
    SMILES string), extracts the coordinate and topology data for the first 
    available structure, and runs the AmberTools `parmchk2` executable. 
    `parmchk2` searches for missing force field parameters (bonds, angles, 
    dihedrals) and appends reasonable parameter suggestions to an output 
    Force Field Modification (`.frcmod`) file.

    Parameters
    ----------
    inp : str, Struct, or ListOfStruct
        Input molecular system. Can be a filepath to a structural file 
        (e.g., Mol2, SMILES), a single `Struct` object, or a `ListOfStruct` 
        collection.
    out : str
        The filepath where the generated `.frcmod` file will be saved.
    clean : bool, default False
        If True, the temporary working directory is deleted upon completion.
        If False, intermediate files are preserved for inspection or debugging.
    verbose : bool, default True
        If True, prints progress updates and temporary workspace directories 
        to standard output.
    parmchk_ff : str, default "gaff"
        The base Amber force field selection used by `parmchk2` to match existing 
        parameters (e.g., 'gaff', 'gaff2'). Corresponds to the `-s` flag.
    parmchk_frc : str, optional
        Filepath to an additional, user-supplied custom `frcmod` parameter file. 
        Corresponds to the `-p` flag.
    parmchk_keep : bool, default False
        If True, tells `parmchk2` to keep intermediate report files (such as 
        missing parameter lists). Corresponds to the `-k` flag.
    parmchk_all : bool, default False
        If True, prints all force field parameters for the molecule into the 
        output file instead of only the missing parameters. Corresponds 
        to the `-a` flag.
    parmchk_atc : str, optional
       This is the name of a supplemental atom type score file. The file is
       used by parmchk to lookup parameters for similar atom types based on
       a score. The standard set of amber and gaff atom types are in a
       global score file ${AMBERHOME}/dat/antechamber/PARMCHK.DAT.
       If you introduce new atom types, you'll either need to modify that
       file or create additional entries in a separate file. The score_file
       is the separate file. It is directly passed to parmchk2 via the -atc
       option. For more information, see "parmchk2 -l"
    parmchk_ffmissing : str, default "gaff"
       The force field parameter set identifier used to lookup missing parameters.
       The only supported values are:
    
        * ``"1"`` or ``"gaff"`` : General Amber Force Field (GAFF)
        * ``"2"`` or ``"gaff2"`` : General Amber Force Field v2 (GAFF2)

        The default is to lookup missing parameters from gaff.
    parmchk_afrc : str, optional
       Additional input frcmod file that supplements the parmchk_ff and parmchk_frc
       parameters.
    only_missing : bool, default False
       If true, then use the frcmod standard of only printing the missing parameters.
       If false (default), then include a MASS and NONBON entry for every type
       appearing in the output frcmod. This allows parmed to read, manipulate, and
       rewrite the frcmod file.
    **standard_kwargs : dict
        Additional keyword arguments evaluated against standard internal model 
        defaults. Unexpected parameters raise a `TypeError`.

    Returns
    -------
    res : dict or output type of run_parmchk2
        The results, return code, or dictionary metadata generated by the 
        underlying `run_parmchk2` wrapper function.

    Raises
    ------
    TypeError
        If an unexpected keyword argument is supplied through `**standard_kwargs`.
    IndexError
        If the processed input generates an empty `ListOfStruct` collection, 
        preventing extraction of the first index (`los[0]`).

    See Also
    --------
    Parmchk.run_parmchk2 : Core Python wrapper executing the `parmchk2` binary.
    Antechamber.run_antechamber : Utility often run prior to generate standard 
                                  Mol2 inputs.

    Examples
    --------
    >>> from my_module import RunCreateAmberFrcmod
    >>> result = RunCreateAmberFrcmod(
    ...     inp="ligand.mol2",
    ...     out="ligand.frcmod",
    ...     parmchk_ff="gaff2",
    ...     parmchk_all=False
    ... )
    """

    import argparse
    from types import SimpleNamespace

    _p = argparse.ArgumentParser(add_help=False)
    #AddModelOptions(_p)
    std_defaults = vars(_p.parse_args([]))
    unknown = set(standard_kwargs) - set(std_defaults)
    if unknown:
        raise TypeError(
            f"Unexpected keyword argument(s): {sorted(unknown)}"
        )
    std = {**std_defaults, **standard_kwargs}

    args = SimpleNamespace(
        inp=inp,
        out=out,
        clean=clean,
        verbose=verbose,
        parmchk_ff=parmchk_ff,
        parmchk_frc=parmchk_frc,
        parmchk_keep=parmchk_keep,
        parmchk_all=parmchk_all,
        parmchk_ffmissing=parmchk_ffmissing,
        parmchk_atc=parmchk_atc,
        parmchk_afrc=parmchk_afrc,
        only_missing=only_missing,
        **std,
    )

    #########################################################################

    from . Parmchk import run_parmchk2
    from . Struct import ListOfStruct, Struct
    import os, copy, tempfile

    if isinstance(args.inp,ListOfStruct):
        los = copy.deepcopy(args.inp)
    elif isinstance(inp,Struct):
        los = ListOfStruct( [args.inp] )
    else:
        los = ListOfStruct.from_file(args.inp)
    
    tmpbase = "tmpdir"
    os.makedirs(tmpbase, exist_ok=True)
    
    # 1. Manually create the unique temporary directory path
    tmpdir_path = tempfile.mkdtemp(dir=tmpbase)

    if args.verbose:
        print(f"Created temporary directory: {tmpdir_path}")
    
    try:

        resname = los[0].data["resnames"][0]
        tmpname = f"{tmpdir_path}/{resname}"
        mol2   = f"{tmpname}.mol2"
        frcmod = args.out
        los[0].SaveCrds( f"{tmpname}.mol2" )

        # Normalize parameters based on internal mapping definitions
        clean_ff = str(parmchk_ff).strip().lower()
        actual_ff = _BASE_FF_MAP.get(clean_ff, {}).get("parmchk", parmchk_ff)
        
        res = run_parmchk2\
            (mol2,frcmod,
             ff_selection=actual_ff,
             ff_frcmod=args.parmchk_frc,
             keep_files=args.parmchk_keep,
             all_parameters=args.parmchk_all,
             verbose=args.verbose,
             ff_missing=args.parmchk_ffmissing,
             score_file=args.parmchk_atc,
             ff_afrc=args.parmchk_afrc,
             only_missing=args.only_missing)

        
    finally:
        # 4. Clean up only if requested
        if clean:
            if args.verbose:
                print(f"\nCleaning up: Removing {tmpdir_path}")
            import shutil
            shutil.rmtree(tmpdir_path)
        elif args.verbose:
            print(f"\n'clean' is False. Leaving directory intact at: {tmpdir_path}")





# def RunCreateAmberTopology(
#         *,
#         mol2_file: str,
#         frcmod_file: str,
#         parmchk_ff: str = "gaff",
#         parmchk_frc: str = None,
#         resname: str = "MOL",
#         out_prefix: str = "system",
#         clean: bool = True,
#         verbose: bool = True
# ):
#     """Run tleap in an isolated temporary directory to construct system topologies.

#     Parameters
#     ----------
#     mol2_file : str, Struct, or ListOfStruct
#         Path to the structural input Mol2 file or an in-memory coordinate structure object.
#     frcmod_file : str
#         Path to the modifications file output by parmchk2.
#     parmchk_ff : str, default "gaff"
#         The identifier passed to parmchk2's `-s` flag. Supports integers ('1'-'5') 
#         or names.
#     parmchk_frc : str, optional
#         Modifier parameters chain split by "+" signs (e.g., 'ff14SB+ol3').
#     resname : str, default "MOL"
#         The residue designation tag embedded within your structural input coordinate files.
#     out_prefix : str, default "system"
#         Naming prefix token applied when writing final library and parameter target outputs.
#     clean : bool, default True
#         If True, the temporary sandbox folder environment is wiped upon completion.
#     verbose : bool, default True
#         If True, outputs active console updates and redirects subprocess streams.

#     Returns
#     -------
#     topology_results : dict
#         Absolute paths mapping the generated assets and the localized script details:
        
#         - "lib" : Absolute path to the saved LEaP object library (.lib/.off) file.
#         - "parm7" : Absolute path to the Amber structural topology (.parm7) file.
#         - "rst7" : Absolute path to the coordinate restart (.rst7) file.
#         - "tleap" : Local text string showing the compiled script assuming files occupy 
#                     the execution directory scope.

#     Raises
#     ------
#     RuntimeError
#         If `tleap` processing crashes or fails to output target system files.

#     See Also
#     --------
#     RunBuildAmberSystem : High-level automation orchestrator for parameters and topologies.
#     """
#     import os
#     import subprocess
#     import shutil
#     import tempfile
#     import traceback
#     import sys
#     from . Struct import ListOfStruct, Struct

#     if verbose:
#         print("[RunCreateAmberTopology] Entering function.")

#     try:
#         dest_dir = os.getcwd()
#         abs_lib_out = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.lib"))
#         abs_parm7_out = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.parm7"))
#         abs_rst7_out = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.rst7"))

#         tmpbase = "tmpdir"
#         os.makedirs(tmpbase, exist_ok=True)
#         tmpdir_path = tempfile.mkdtemp(dir=tmpbase, prefix="tleap_")
        
#         if verbose:
#             print(f"[RunCreateAmberTopology] Created isolated workspace: {tmpdir_path}")

#         local_mol2 = f"{out_prefix}_topology_input.mol2"
#         local_frcmod = os.path.basename(frcmod_file)
#         local_lib = f"{out_prefix}.lib"
#         local_parm7 = f"{out_prefix}.parm7"
#         local_rst7 = f"{out_prefix}.rst7"
#         script_name = "tleap_commands.in"

#         if isinstance(mol2_file, (Struct, ListOfStruct)):
#             target_mol2_path = os.path.join(tmpdir_path, local_mol2)
#             if verbose:
#                 print(f"[RunCreateAmberTopology] Exporting in-memory object data via SaveCrds: {target_mol2_path}")
#             mol2_file.SaveCrds(target_mol2_path)
#         else:
#             local_mol2 = os.path.basename(str(mol2_file))
#             shutil.copy2(str(mol2_file), os.path.join(tmpdir_path, local_mol2))

#         shutil.copy2(frcmod_file, os.path.join(tmpdir_path, local_frcmod))

#         tleap_commands = ["# Automatically mapped options framework"]
#         tleap_commands.extend(resolve_tleap_commands(parmchk_ff, parmchk_frc))
#         tleap_commands.append(f"loadAmberParams {local_frcmod}")
#         tleap_commands.extend([
#             f"{resname} = loadMol2 {local_mol2}",
#             f"saveOff {resname} {local_lib}",
#             f"saveAmberParm {resname} {local_parm7} {local_rst7}",
#             "quit"
#         ])

#         tleap_script_content = "\n".join(tleap_commands) + "\n"

#         with open(os.path.join(tmpdir_path, script_name), "w") as f:
#             f.write(tleap_script_content)

#         try:
#             cmd = ["tleap", "-f", script_name]
#             result = subprocess.run(
#                 cmd, 
#                 stdout=subprocess.PIPE, 
#                 stderr=subprocess.PIPE, 
#                 text=True, 
#                 cwd=tmpdir_path
#             )
            
#             if verbose:
#                 print(result.stdout)
#                 if result.stderr:
#                     print(f"[RunCreateAmberTopology] Subprocess Stderr Output:\n{result.stderr}")
                
#             tmp_parm7 = os.path.join(tmpdir_path, local_parm7)
#             tmp_rst7 = os.path.join(tmpdir_path, local_rst7)
#             tmp_lib = os.path.join(tmpdir_path, local_lib)

#             if not os.path.exists(tmp_parm7):
#                 err_msg = (
#                     f"tleap failed to generate the parameter topology file '{local_parm7}'. "
#                     f"This typically indicates missing atom types in your force field configurations. "
#                     f"Subprocess exit code: {result.returncode}."
#                 )
#                 raise RuntimeError(err_msg)
                
#             shutil.move(tmp_lib, abs_lib_out)
#             shutil.move(tmp_parm7, abs_parm7_out)
#             shutil.move(tmp_rst7, abs_rst7_out)

#             if verbose:
#                 print("[RunCreateAmberTopology] Exiting function: Completed successfully.")

#             return {
#                 "lib": abs_lib_out,
#                 "parm7": abs_parm7_out,
#                 "rst7": abs_rst7_out,
#                 "tleap": tleap_script_content
#             }
            
#         finally:
#             if clean:
#                 if verbose:
#                     print(f"[RunCreateAmberTopology] Scrubbing sandbox directory: {tmpdir_path}")
#                 shutil.rmtree(tmpdir_path)

#     except Exception as e:
#         print(f"[RunCreateAmberTopology] Error: Function execution failed. Reason: {str(e)}", file=sys.stderr)
#         print("-" * 60, file=sys.stderr)
#         traceback.print_exc(file=sys.stderr)
#         print("-" * 60, file=sys.stderr)
#         raise


# def RunBuildAmberSystem(
#         *,
#         inp: str,
#         parmchk_ff: str = "gaff",
#         parmchk_frc: str = None,
#         clean: bool = True,
#         verbose: bool = True,
#         tar: bool = False,
#         **standard_kwargs
# ):
#     """Orchestrate parmchk2 and tleap pipelines using an input molecule configuration.

#     Parameters
#     ----------
#     inp : str or Struct
#         Molecular geometry input. Can be a system file path string (e.g., Mol2) 
#         or an active in-memory `Struct` instance.
#     parmchk_ff : str, default "gaff"
#         Shared configuration baseline force field selection token passed directly to parmchk2.
#     parmchk_frc : str, optional
#         Additional parameter modifier tokens passed directly to parmchk2's `-frc` parameter.
#     clean : bool, default True
#         If True, deletes individual transient sandbox folders and auxiliary files.
#     verbose : bool, default True
#         If True, traces pipeline execution milestones to standard output.
#     tar : bool, default False
#         If True, compresses all resulting outputs safely into an isolated 
#         '{resname}.tar.gz' file alongside a reproducible reconstruction shell script.
#     **standard_kwargs : dict
#         Additional keyword arguments evaluated against standard internal model defaults.

#     Returns
#     -------
#     topology_results : dict
#         Absolute local system file paths mapping out the generated files, the target 
#         script text, and optionally the generated archive:
        
#         - "lib" : Absolute path to the saved LEaP object library (.lib/.off) file.
#         - "parm7" : Absolute path to the Amber structural topology (.parm7) file.
#         - "rst7" : Absolute path to the coordinate restart (.rst7) file.
#         - "tleap" : Localized text string representing the executed tleap script.
#         - "tar" : Absolute path to the archive file, included only if `tar=True`.

#     Raises
#     ------
#     TypeError
#         If an unexpected keyword argument is supplied through `**standard_kwargs`.
#     RuntimeError
#         If parmchk2 or tleap crashes or catches missing atom type parameters.

#     See Also
#     --------
#     RunCreateAmberFrcmod : Generates the missing parameter modification file.
#     RunCreateAmberTopology : Handles the isolated execution layer of tleap.
#     """
#     import os
#     import shutil
#     import tarfile
#     import tempfile
#     import traceback
#     import sys
#     from . Struct import ListOfStruct, Struct

#     if verbose:
#         print("[RunBuildAmberSystem] Entering function.")

#     try:
#         # 1. Resolve the Struct instance 's' directly from the variant 'inp' type definitions
#         if isinstance(inp, Struct):
#             if verbose:
#                 print("[RunBuildAmberSystem] Input is a Struct object. Creating a local copy.")
#             s = inp.copy()
#         else:
#             if verbose:
#                 print(f"[RunBuildAmberSystem] Input is a filepath string: {inp}. Loading ListOfStruct.")
#             los = ListOfStruct.from_file(str(inp))
#             s = los

#         # 2. Extract out_prefix and resname parameters from s.data["resnames"][0]
#         resname = str(s.data["resnames"][0]).strip()
#         out_prefix = resname
        
#         if verbose:
#             print(f"[RunBuildAmberSystem] Extracted active resname: '{resname}' -> out_prefix: '{out_prefix}'")

#         # Establish top-level thread-safe root scratch folder tracking 
#         tmpbase = "tmpdir"
#         os.makedirs(tmpbase, exist_ok=True)
#         orchestrator_scratch = tempfile.mkdtemp(dir=tmpbase, prefix=f"build_{resname}_")

#         if verbose:
#             print(f"[RunBuildAmberSystem] Created thread-isolated directory: {orchestrator_scratch}")

#         # Build unique internal relative paths using the clean f"{out_prefix}.mol2" layout
#         local_frcmod = os.path.join(orchestrator_scratch, f"{out_prefix}.frcmod")
#         mol2_source = os.path.join(orchestrator_scratch, f"{out_prefix}.mol2")
        
#         if verbose:
#             print(f"[RunBuildAmberSystem] Exporting coordinates data inside scratchpad via SaveCrds: {mol2_source}")
#         s.SaveCrds(mol2_source)

#         if verbose:
#             print("[RunBuildAmberSystem] >>> Step 1: Identifying parameters dependencies via parmchk2...")
            
#         RunCreateAmberFrcmod(
#             inp=mol2_source,
#             out=local_frcmod,
#             clean=clean,
#             verbose=verbose,
#             parmchk_ff=parmchk_ff,
#             parmchk_frc=parmchk_frc,
#             **standard_kwargs
#         )
        
#         if verbose:
#             print("[RunBuildAmberSystem] >>> Step 2: Compiling topological structures via tleap...")

#         # RunCreateAmberTopology natively saves out the uncompressed .lib, .parm7, and .rst7 files to the active directory
#         topology_results = RunCreateAmberTopology(
#             mol2_file=mol2_source,
#             frcmod_file=local_frcmod,
#             parmchk_ff=parmchk_ff,
#             parmchk_frc=parmchk_frc,
#             resname=resname,
#             out_prefix=out_prefix,
#             clean=clean,
#             verbose=verbose
#         )
        
#         # Track our destination workspace directory
#         dest_dir = os.getcwd()
        
#         # Ensure that the loose .mol2 and .frcmod files are also copied into the current directory
#         abs_mol2_dest = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.mol2"))
#         abs_frcmod_dest = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.frcmod"))
#         shutil.copy2(mol2_source, abs_mol2_dest)
#         shutil.copy2(local_frcmod, abs_frcmod_dest)
        
#         # Write out the loose version of the executed tleap commands file locally
#         final_tleap_in_path = os.path.abspath(os.path.join(dest_dir, "tleap_commands.in"))
#         local_tleap_in = os.path.join(orchestrator_scratch, "tleap_commands.in")
#         with open(local_tleap_in, "w") as sf:
#             sf.write(topology_results["tleap"])
#         shutil.copy2(local_tleap_in, final_tleap_in_path)

#         # Generate the executable reproduction shell script text
#         bash_mol2_name = os.path.basename(mol2_source)
#         bash_frcmod_name = f"{out_prefix}.frcmod"
        
#         clean_ff = str(parmchk_ff).strip().lower()
#         actual_ff = _BASE_FF_MAP.get(clean_ff, {}).get("parmchk", parmchk_ff)
        
#         parmchk_cmd = f"parmchk2 -i {bash_mol2_name} -f mol2 -o {bash_frcmod_name} -s {actual_ff}"
#         if parmchk_frc:
#             parmchk_cmd += f" -frc {parmchk_frc}"
            
#         bash_script_content = f"""#!/bin/bash
# # Automatically generated reproduction execution pipeline script
# # Reconstructs topology and parameters assets from raw structural coordinates

# echo ">>> Re-generating modification parameter records via parmchk2..."
# {parmchk_cmd}

# echo ">>> Re-compiling system simulation objects via tleap..."
# tleap -f tleap_commands.in

# echo ">>> Structural reconstruction completed successfully."
# """
        
#         local_bash_sh = os.path.join(orchestrator_scratch, f"{out_prefix}.sh")
#         with open(local_bash_sh, "w") as bf:
#             bf.write(bash_script_content)
#         os.chmod(local_bash_sh, 0o755)
        
#         # Ensure that the loose reconstruction script is also placed in the current directory
#         final_bash_sh_path = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.sh"))
#         shutil.copy2(local_bash_sh, final_bash_sh_path)
            
#         if tar:
#             final_tar_path = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.tar.gz"))
#             try:
#                 if verbose:
#                     print(f"[RunBuildAmberSystem] >>> Step 3: Packing generated system assets into archive {final_tar_path}...")
                    
#                 # Compile archive inside scratchpad using local paths before publishing out to workspace root
#                 temp_tar_path = os.path.join(orchestrator_scratch, f"{out_prefix}.tar.gz")
#                 with tarfile.open(temp_tar_path, "w:gz") as tar_file:
#                     tar_file.add(mol2_source, arcname=bash_mol2_name)
#                     tar_file.add(local_frcmod, arcname=bash_frcmod_name)
#                     tar_file.add(topology_results["lib"], arcname=os.path.basename(topology_results["lib"]))
#                     tar_file.add(topology_results["parm7"], arcname=os.path.basename(topology_results["parm7"]))
#                     tar_file.add(topology_results["rst7"], arcname=os.path.basename(topology_results["rst7"]))
#                     tar_file.add(local_tleap_in, arcname="tleap_commands.in")
#                     tar_file.add(local_bash_sh, arcname=f"{out_prefix}.sh")
                    
#                 # Publish the built compressed package cleanly to active working space root
#                 shutil.move(temp_tar_path, final_tar_path)
#                 topology_results["tar"] = final_tar_path
                
#             except Exception as tar_err:
#                 raise RuntimeError(f"Failed to generate or write out file assets to tar compression archive: {str(tar_err)}")

#         if clean:
#             if verbose:
#                 print(f"[RunBuildAmberSystem] Removing scratchpad workspace tracking directory: {orchestrator_scratch}")
#             shutil.rmtree(orchestrator_scratch)
            
#         if verbose:
#             print("[RunBuildAmberSystem] Exiting function: Completed successfully.")
            
#         return topology_results

#     except Exception as e:
#         print(f"[RunBuildAmberSystem] Error: Orchestration pipeline wrapper failed. Reason: {str(e)}", file=sys.stderr)
#         print("-" * 60, file=sys.stderr)
#         traceback.print_exc(file=sys.stderr)
#         print("-" * 60, file=sys.stderr)
#         raise













def RunCreateAmberTopology(
        *,
        mol2_file: str,
        frcmod_file: str = None,
        parmchk_ff: str = "gaff",
        parmchk_frc: str = None,
        resname: str = "MOL",
        out_prefix: str = "system",
        clean: bool = True,
        verbose: bool = True
):
    """Run tleap in an isolated temporary directory to construct system topologies.

    Parameters
    ----------
    mol2_file : str, Struct, or ListOfStruct
        Path to the structural input Mol2 file or an in-memory coordinate structure object.
    frcmod_file : str or list of str, optional, default None
        Path string or list of path strings pointing to force field modification (.frcmod) files.
        If None, no modification files are loaded into the simulation workspace macro.
    parmchk_ff : str, default "gaff"
        The identifier passed to parmchk2's `-s` flag. Supports integers ('1'-'5') or names.
    parmchk_frc : str, optional
        Modifier parameters chain split by "+" signs (e.g., 'ff14SB+ol3').
    resname : str, default "MOL"
        The residue designation tag embedded within your structural input coordinate files.
    out_prefix : str, default "system"
        Naming prefix token applied when writing final library and parameter target outputs.
    clean : bool, default True
        If True, the temporary sandbox folder environment is wiped upon completion.
    verbose : bool, default True
        If True, outputs active console updates and redirects subprocess streams.

    Returns
    -------
    topology_results : dict
        Absolute paths mapping the generated assets and the localized script details:
        
        - "lib" : Absolute path to the saved LEaP object library (.lib/.off) file.
        - "parm7" : Absolute path to the Amber structural topology (.parm7) file.
        - "rst7" : Absolute path to the coordinate restart (.rst7) file.
        - "tleap" : Local text string showing the compiled script assuming files occupy 
                    the execution directory scope.

    Raises
    ------
    RuntimeError
        If `tleap` processing crashes or fails to output target system files.
    """
    import os
    import subprocess
    import shutil
    import tempfile
    import traceback
    import sys
    from . Struct import ListOfStruct, Struct

    if verbose:
        print("[RunCreateAmberTopology] Entering function.")

    try:
        # Standardize our input variants (None, str, list) into a uniform array list loop
        if frcmod_file is None:
            frcmod_list = []
        elif isinstance(frcmod_file, str):
            frcmod_list = [frcmod_file]
        elif isinstance(frcmod_file, list):
            frcmod_list = frcmod_file
        else:
            frcmod_list = list(frcmod_file) if frcmod_file else []

        dest_dir = os.getcwd()
        abs_lib_out = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.lib"))
        abs_parm7_out = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.parm7"))
        abs_rst7_out = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.rst7"))

        tmpbase = "tmpdir"
        os.makedirs(tmpbase, exist_ok=True)
        tmpdir_path = tempfile.mkdtemp(dir=tmpbase, prefix="tleap_")
        
        if verbose:
            print(f"[RunCreateAmberTopology] Created isolated workspace: {tmpdir_path}")

        local_mol2 = f"{out_prefix}_topology_input.mol2"
        local_lib = f"{out_prefix}.lib"
        local_parm7 = f"{out_prefix}.parm7"
        local_rst7 = f"{out_prefix}.rst7"
        script_name = "tleap_commands.in"

        if isinstance(mol2_file, (Struct, ListOfStruct)):
            target_mol2_path = os.path.join(tmpdir_path, local_mol2)
            if verbose:
                print(f"[RunCreateAmberTopology] Exporting in-memory object data via SaveCrds: {target_mol2_path}")
            mol2_file.SaveCrds(target_mol2_path)
        else:
            local_mol2 = os.path.basename(str(mol2_file))
            shutil.copy2(str(mol2_file), os.path.join(tmpdir_path, local_mol2))

        # Build commands array referencing the globally scoped helper
        tleap_commands = ["# Automatically mapped options framework"]
        tleap_commands.extend(resolve_tleap_commands(parmchk_ff, parmchk_frc))
        
        # Iteratively copy and append sequential load instructions for every resolved modification target
        for f_path in frcmod_list:
            local_frcmod_name = os.path.basename(f_path)
            shutil.copy2(f_path, os.path.join(tmpdir_path, local_frcmod_name))
            tleap_commands.append(f"loadAmberParams {local_frcmod_name}")

        tleap_commands.extend([
            f"{resname} = loadMol2 {local_mol2}",
            f"saveOff {resname} {local_lib}",
            f"saveAmberParm {resname} {local_parm7} {local_rst7}",
            "quit"
        ])

        tleap_script_content = "\n".join(tleap_commands) + "\n"

        with open(os.path.join(tmpdir_path, script_name), "w") as f:
            f.write(tleap_script_content)

        try:
            cmd = ["tleap", "-f", script_name]
            
            if verbose:
                print("[RunCreateAmberTopology] tleap script %s"%(os.path.join(tmpdir_path, script_name)))
                print(tleap_script_content)
                print("[RunCreateAmberTopology] Running: %s"%(" ".join(cmd)))
            
            result = subprocess.run(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                text=True, 
                cwd=tmpdir_path
            )
            
            if verbose:
                print(result.stdout)
                if result.stderr:
                    print(f"[RunCreateAmberTopology] Subprocess Stderr Output:\n{result.stderr}")
                
            tmp_parm7 = os.path.join(tmpdir_path, local_parm7)
            tmp_rst7 = os.path.join(tmpdir_path, local_rst7)
            tmp_lib = os.path.join(tmpdir_path, local_lib)

            if not os.path.exists(tmp_parm7):
                err_msg = (
                    f"tleap failed to generate the parameter topology file '{local_parm7}'. "
                    f"Subprocess exit code: {result.returncode}."
                )
                raise RuntimeError(err_msg)
                
            shutil.move(tmp_lib, abs_lib_out)
            shutil.move(tmp_parm7, abs_parm7_out)
            shutil.move(tmp_rst7, abs_rst7_out)

            if verbose:
                print("[RunCreateAmberTopology] Exiting function: Completed successfully.")

            return {
                "lib": abs_lib_out,
                "parm7": abs_parm7_out,
                "rst7": abs_rst7_out,
                "tleap": tleap_script_content
            }
            
        finally:
            if clean:
                if verbose:
                    print(f"[RunCreateAmberTopology] Scrubbing sandbox directory: {tmpdir_path}")
                shutil.rmtree(tmpdir_path)

    except Exception as e:
        print(f"[RunCreateAmberTopology] Error: Function execution failed. Reason: {str(e)}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        raise

    







def RunBuildAmberSystem(
        *,
        inp: str,
        parmchk_ff: str = "gaff",
        parmchk_frc: str = None,
        frcmod_files: list = None,
        clean: bool = True,
        verbose: bool = True,
        tar: bool = False,
        parmchk_ffmissing: str = "gaff",
        **standard_kwargs
):
    """Orchestrate parmchk2 and tleap pipelines using an input molecule configuration.

    Parameters
    ----------
    inp : str or Struct
        Molecular geometry input. Can be a system file path string (e.g., Mol2) 
        or an active in-memory `Struct` instance.
    parmchk_ff : str, default "gaff"
        Shared configuration baseline force field selection token passed directly to parmchk2.
    parmchk_frc : str, optional
        Additional parameter modifier tokens passed directly to parmchk2's `-frc` parameter.
    frcmod_files : str or list of str, optional, default None
        Explicit paths to pre-generated .frcmod files. Can be a single string path or a list. 
        If None, the pipeline dynamically creates a missing parameter file via parmchk2.
    clean : bool, default True
        If True, deletes individual transient sandbox folders and auxiliary files.
    verbose : bool, default True
        If True, traces pipeline execution milestones to standard output.
    tar : bool, default False
        If True, compresses all resulting outputs safely into an isolated 
        '{resname}.tar.gz' file alongside a reproducible reconstruction shell script.
    parmchk_ffmissing : str, default "gaff"
       The force field parameter set identifier used to lookup missing parameters.
       The only supported values are:
    
        * ``"1"`` or ``"gaff"`` : General Amber Force Field (GAFF)
        * ``"2"`` or ``"gaff2"`` : General Amber Force Field v2 (GAFF2)
    
        The default is to lookup missing parameters from gaff.
    **standard_kwargs : dict
        Additional keyword arguments evaluated against standard internal model defaults.

    Returns
    -------
    topology_results : dict
        Absolute local system file paths mapping out the generated files, the target 
        script text, and optionally the generated archive.
    """
    import os
    import shutil
    import tarfile
    import tempfile
    import traceback
    import sys
    from . Struct import ListOfStruct, Struct

    if verbose:
        print("[RunBuildAmberSystem] Entering function.")

    try:
        # Standardize the input parameter variations (None, single str path, or lists) 
        if frcmod_files is None:
            injected_frcmod_list = []
        elif isinstance(frcmod_files, str):
            injected_frcmod_list = [frcmod_files]
        elif isinstance(frcmod_files, list):
            injected_frcmod_list = frcmod_files
        else:
            injected_frcmod_list = list(frcmod_files) if frcmod_files else []

        # 1. Resolve the Struct instance 's' directly from the variant 'inp' type definitions
        if isinstance(inp, Struct):
            if verbose:
                print("[RunBuildAmberSystem] Input is a Struct object. Creating a local copy.")
            s = inp.copy()
        else:
            if verbose:
                print(f"[RunBuildAmberSystem] Input is a filepath string: {inp}. Loading ListOfStruct.")
            los = ListOfStruct.from_file(str(inp))
            s = los[0]

        # 2. Extract out_prefix and resname parameters from s.data["resnames"]
        resname = str(s.data["resnames"][0]).strip()
        out_prefix = resname
        
        if verbose:
            print(f"[RunBuildAmberSystem] Extracted active resname: '{resname}' -> out_prefix: '{out_prefix}'")

        # Establish top-level thread-safe root scratch folder tracking 
        tmpbase = "tmpdir"
        os.makedirs(tmpbase, exist_ok=True)
        orchestrator_scratch = tempfile.mkdtemp(dir=tmpbase, prefix=f"build_{resname}_")

        if verbose:
            print(f"[RunBuildAmberSystem] Created thread-isolated directory: {orchestrator_scratch}")

        # Build unique internal paths inside the local thread scratchpad space
        mol2_source = os.path.join(orchestrator_scratch, f"{out_prefix}.mol2")
        if verbose:
            print(f"[RunBuildAmberSystem] Exporting coordinates data inside scratchpad via SaveCrds: {mol2_source}")
        s.SaveCrds(mol2_source)

        dest_dir = os.getcwd()
        resolved_frcmod_targets = []

        # 3. Determine if we skip parameter discovery or load explicit files from our list
        if len(injected_frcmod_list) > 0:
            if verbose:
                print(f"[RunBuildAmberSystem] Explicit frcmod file(s) input detected: {injected_frcmod_list}. Skipping parameter creation.")
            
            for user_frcmod in injected_frcmod_list:
                if not os.path.exists(user_frcmod):
                    raise FileNotFoundError(f"Injected user frcmod file path not found: {user_frcmod}")
                
                scratch_frcmod_copy = os.path.join(orchestrator_scratch, os.path.basename(user_frcmod))
                shutil.copy2(user_frcmod, scratch_frcmod_copy)
                resolved_frcmod_targets.append(scratch_frcmod_copy)
                
                # Mirror them out to the current working public workspace directory immediately
                public_frcmod_mirror = os.path.abspath(os.path.join(dest_dir, os.path.basename(user_frcmod)))
                src = Path(user_frcmod)
                dst = Path(public_frcmod_mirror)
                #shutil.copy2(user_frcmod, public_frcmod_mirror)
                if src.exists() and dst.exists() and src.samefile(dst):
                    pass  # Skip copying because they are identical
                else:
                    shutil.copy2(src, dst)
        else:
            if verbose:
                print("[RunBuildAmberSystem] >>> Step 1: Identifying parameters dependencies via parmchk2...")
            
            local_frcmod = os.path.join(orchestrator_scratch, f"{out_prefix}.frcmod")
            RunCreateAmberFrcmod(
                inp=mol2_source,
                out=local_frcmod,
                clean=clean,
                verbose=verbose,
                parmchk_ff=parmchk_ff,
                parmchk_frc=parmchk_frc,
                parmchk_ffmissing=parmchk_ffmissing,
                **standard_kwargs
            )
            resolved_frcmod_targets.append(local_frcmod)
            
            # Mirror the parmchk2 generated frcmod out to the current working public directory
            abs_frcmod_dest = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.frcmod"))
            shutil.copy2(local_frcmod, abs_frcmod_dest)

        if verbose:
            print("[RunBuildAmberSystem] >>> Step 2: Compiling topological structures via tleap...")

        # Pass our normalized array collection downstream directly into our list-aware compiler
        topology_results = RunCreateAmberTopology(
            mol2_file=mol2_source,
            frcmod_file=resolved_frcmod_targets if resolved_frcmod_targets else None,
            parmchk_ff=parmchk_ff,
            parmchk_frc=parmchk_frc,
            resname=resname,
            out_prefix=out_prefix,
            clean=clean,
            verbose=verbose
        )
        
        # Ensure that the loose .mol2, reproduction script, and tleap configurations match exactly
        abs_mol2_dest = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.mol2"))
        shutil.copy2(mol2_source, abs_mol2_dest)
        
        final_tleap_in_path = os.path.abspath(os.path.join(dest_dir, "tleap_commands.in"))
        local_tleap_in = os.path.join(orchestrator_scratch, "tleap_commands.in")
        with open(local_tleap_in, "w") as sf:
            sf.write(topology_results["tleap"])
        shutil.copy2(local_tleap_in, final_tleap_in_path)

        bash_mol2_name = os.path.basename(mol2_source)
        
        if len(injected_frcmod_list) > 0:
            parmchk_cmd = "# Skipping parmchk2 step: loaded explicit pre-existing frcmod targets file list"
        else:
            bash_frcmod_name = f"{out_prefix}.frcmod"

            ############################################################
            # Map parameters to ffpopt-parmchk2.py CLI flags
            cmd_parts = [
                f"ffpopt-parmchk2.py",
                f"-i {bash_mol2_name}",
                f"-o {bash_frcmod_name}",
                f"-s {parmchk_ff}"
            ]
            
            if parmchk_frc:
                cmd_parts.append(f"-frc {parmchk_frc}")
            if parmchk_ffmissing:
                cmd_parts.append(f"--ffmiss {parmchk_ffmissing}")
            if verbose:
                cmd_parts.append("--verbose")
            if clean:
                cmd_parts.append("--clean")
                
            parmchk_cmd = " ".join(cmd_parts)

        bash_script_content = f"""#!/bin/bash
# Automatically generated reproduction execution pipeline script
# Reconstructs topology and parameters assets from raw structural coordinates

echo ">>> Setting up parameter traces..."
{parmchk_cmd}

echo ">>> Re-compiling system simulation objects via tleap..."
tleap -f tleap_commands.in

echo ">>> Structural reconstruction completed successfully."
"""
        ############################################################
       
#             clean_ff = str(parmchk_ff).strip().lower()
#             actual_ff = _BASE_FF_MAP.get(clean_ff, {}).get("parmchk", parmchk_ff)
#             parmchk_cmd = f"parmchk2 -i {bash_mol2_name} -f mol2 -o {bash_frcmod_name} -s {actual_ff}"
#             if parmchk_frc:
#                 parmchk_cmd += f" -frc {parmchk_frc}"
            
#         bash_script_content = f"""#!/bin/bash
# # Automatically generated reproduction execution pipeline script
# # Reconstructs topology and parameters assets from raw structural coordinates

# echo ">>> Setting up parameter traces..."
# {parmchk_cmd}

# echo ">>> Re-compiling system simulation objects via tleap..."
# tleap -f tleap_commands.in

# echo ">>> Structural reconstruction completed successfully."
# """
        
        local_bash_sh = os.path.join(orchestrator_scratch, f"{out_prefix}.sh")
        with open(local_bash_sh, "w") as bf:
            bf.write(bash_script_content)
        os.chmod(local_bash_sh, 0o755)
        
        final_bash_sh_path = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.sh"))
        shutil.copy2(local_bash_sh, final_bash_sh_path)
            
        if tar:
            final_tar_path = os.path.abspath(os.path.join(dest_dir, f"{out_prefix}.tar.gz"))
            try:
                if verbose:
                    print(f"[RunBuildAmberSystem] >>> Step 3: Packing generated system assets into archive {final_tar_path}...")
                    
                temp_tar_path = os.path.join(orchestrator_scratch, f"{out_prefix}.tar.gz")
                with tarfile.open(temp_tar_path, "w:gz") as tar_file:
                    tar_file.add(mol2_source, arcname=bash_mol2_name)
                    tar_file.add(topology_results["lib"], arcname=os.path.basename(topology_results["lib"]))
                    tar_file.add(topology_results["parm7"], arcname=os.path.basename(topology_results["parm7"]))
                    tar_file.add(topology_results["rst7"], arcname=os.path.basename(topology_results["rst7"]))
                    tar_file.add(local_tleap_in, arcname="tleap_commands.in")
                    tar_file.add(local_bash_sh, arcname=f"{out_prefix}.sh")
                    
                    # Package each explicit frcmod file asset cleanly into the tarball root
                    for scratch_frcmod in resolved_frcmod_targets:
                        tar_file.add(scratch_frcmod, arcname=os.path.basename(scratch_frcmod))
                    
                shutil.move(temp_tar_path, final_tar_path)
                topology_results["tar"] = final_tar_path
            except Exception as tar_err:
                raise RuntimeError(f"Failed to generate or write out file assets to tar compression archive: {str(tar_err)}")

        if clean:
            if verbose:
                print(f"[RunBuildAmberSystem] Removing scratchpad workspace tracking directory: {orchestrator_scratch}")
            shutil.rmtree(orchestrator_scratch)
            
        if verbose:
            print("[RunBuildAmberSystem] Exiting function: Completed successfully.")
            
        return topology_results

    except Exception as e:
        print(f"[RunBuildAmberSystem] Error: Orchestration pipeline wrapper failed. Reason: {str(e)}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        raise



