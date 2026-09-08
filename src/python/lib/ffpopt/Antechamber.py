"""
Antechamber Consolidated Interface Module.

This module provides a unified interface to AmberTools' antechamber program to
correct atom names, assign force field atom types, and calculate partial charges
within a controlled filesystem sandbox using ParmEd and OpenBabel.

Notes
-----
Input structural preparations utilize OpenBabel pybel to perceive bond topology 
graphs in RAM, writing the results directly to the required disk scratch file
to satisfy ParmEd's string file-path requirements natively.
"""

from typing import List, Union, Dict, Optional, Tuple, Any
import numpy as np

def _check_antechamber_executable() -> str:
    """
    Verify the presence of antechamber in system path environments.

    Returns
    -------
    exe_path : str
        The localized file track targeting the functional binary.

    Raises
    ------
    RuntimeError
        If antechamber cannot be located or executed across system paths.
    """
    import shutil
    
    exe = shutil.which("antechamber")
    if exe is None:
        raise RuntimeError(
            "Critical Error: 'antechamber' execution binary was not found "
            "within active environment variables or system path configuration."
        )
    return exe


def FixCharges(inpcharges: List[float], digits: int = 4, target: Optional[int] = None) -> List[float]:
    """
    Adjust and round atomic partial charges to sum precisely to the target or nearest integer.

    Parameters
    ----------
    inpcharges : list of float
        The raw floating-point atomic partial charges to be corrected.
    digits : int, default 4
        The number of decimal places to round the individual values.
    target : int, default None
        Optional designated integer net charge. If None, it is calculated from the inputs.

    Returns
    -------
    fixed_charges : list of float
        The normalized partial atomic charges rounded and balanced to the integer net sum.
    """
    from collections import defaultdict as ddict

    charges = [q for q in inpcharges]
    q = sum(charges)
    nat = len(charges)
    
    if target is None:
        intq = int(round(q))
    else:
        intq = target
        
    dq = (intq - q) / nat
    
    for i in range(nat):
        charges[i] += dq
        charges[i] = float(f"%.{digits}f" % (charges[i]))
        
    qmap = ddict(list)
    for i in range(nat):
        qmap[charges[i]].append(i)

    uniqueqs = [qk for qk in qmap]
    degens = [len(qmap[qk]) for qk in qmap]

    q = sum([u * g for u, g in zip(uniqueqs, degens)])
    sorted_degens = sorted(degens, reverse=True)

    tmpqs = [charge for charge in charges]
    for g in sorted_degens:
        tmpqs = [charge for charge in charges]
        for i in range(len(uniqueqs)):
            if degens[i] == g:
                for j in qmap[uniqueqs[i]]:
                    dq_val = float(f"%.{digits}f" % ((q - intq) / g))
                    tmpqs[j] -= dq_val
                break
        if abs(sum(tmpqs) - intq) < 1e-6 or sum(tmpqs) == 0:
            break
        
    for i in range(nat):
        charges[i] = tmpqs[i]
        
    return charges


def _generate_high_precision_input_mol2(
    coords: Union[List[List[float]], np.ndarray], 
    elements: List[str],
    base_dir: str,
    charge: float,
    log_buffer: Any,
    resname: Optional[str] = "MOL"
) -> str:
    """
    Build a high-precision input MOL2 track file using in-memory pybel routines.

    Parameters
    ----------
    coords : list of list of float or numpy.ndarray
        Shape (N, 3) matrix mapping precise Cartesian coordinates in Angstroms.
    elements : list of str
        Ordered alphanumeric atomic element symbols matching coordinate rows.
    base_dir : str
        Target root folder where the physical coordinate scratch file is saved.
    charge : float
        The absolute formal net charge used to evenly pre-charge the layout atoms.
    log_buffer : io.StringIO
        Memory character buffer used to capture operational diagnostic traces.
    resname : list of str, default "MOL"
        Name of the residue

    Returns
    -------
    temp_path : str
        The secure file path targeting the newly created scratch format.
    """
    import os
    import tempfile
    import parmed as pmd
    from openbabel import pybel
    from collections import defaultdict as ddict
    
    if not os.path.exists(base_dir):
        os.makedirs(base_dir, exist_ok=True)
        
    fd, path = tempfile.mkstemp(suffix=".mol2", dir=base_dir)
    os.close(fd)
    
    num_atoms = len(elements)
    xyz_lines = [f"{num_atoms}", "Generated inside memory via string concatenation"]
    for i in range(num_atoms):
        elem = elements[i]
        x, y, z = coords[i]
        xyz_lines.append(f"{elem:<2} {x:15.8f} {y:15.8f} {z:15.8f}")
    raw_xyz_text = "\n".join(xyz_lines)
    
    pybel_mol = pybel.readstring("xyz", raw_xyz_text)
    
    log_buffer.write("\n======================================================================\n")
    log_buffer.write("HYPOTHESIS TESTING: FILE-BASED PERCEPTION CONVERSION DIAGNOSTICS\n")
    log_buffer.write("======================================================================\n")
    log_buffer.write(f"Total element nodes registered in Pybel block: {len(pybel_mol.atoms)}\n")
    
    try:
        perceived_bonds = pybel_mol.OBMol.NumBonds()
        log_buffer.write(f"OpenBabel perceived bond graph count from layout: {perceived_bonds}\n")
    except Exception as err:
        log_buffer.write(f"OBMol.NumBonds check failed: {err}\n")
    log_buffer.write(f"Writing parsed molecular data straight to target path string: {path}\n")
    log_buffer.write("======================================================================\n\n")

    pybel_mol.write("mol2", path, overwrite=True)
    pmd_struct = pmd.load_file(path, structure=True)

    dq = float(charge) / len(pmd_struct.atoms)
    for atom in pmd_struct.atoms:
        atom.charge = dq
        atom.residue.name = resname
        
    # has_residues = (residxs is not None) and (resnames is not None)
    # if has_residues:
    #     from parmed.topologyobjects import Residue
    #     residue_cache = {}
    #     for i, atom in enumerate(pmd_struct.atoms):
    #         r_idx = residxs[i]
    #         r_name = resnames[i]
    #         if r_idx not in residue_cache:
    #             residue_cache[r_idx] = Residue(name=r_name, number=r_idx + 1)
    #         target_res = residue_cache[r_idx]
            
    #         atom.residue = target_res
    #         target_res.add_atom(atom)
    #         if target_res not in pmd_struct.residues:
    #             pmd_struct.residues.append(target_res)
                
    seen = ddict(int)
    for atom in pmd_struct.atoms:
        atom.type = atom.name
        seen[atom.name] += 1
        count = seen[atom.name]
        if count > 1:
            atom.name = f"{atom.name}{count}"
            
    # from pathlib import Path
    # fpath = Path(path)
    # fpath.unlink(missing_ok=True)
    

    pmd_struct.save(path, format="MOL2", overwrite=True)

    return path


def run_antechamber(
    coords: Union[List[List[float]], np.ndarray],
    elements: List[str],
    charge: int = 0,
    get_charges: bool = True,
    get_types: bool = True,
    get_names: bool = True,
    charge_method: str = "bcc",
    atom_type_family: str = "gaff",
    optimize_geometry: bool = False,
    adjust_names: Optional[str] = None,
    base_dir: str = "tmpfiles",
    resname: Optional[str] = "MOL",
    verbose: bool = False
) -> Dict[str, Optional[List[Union[str, float]]]]:
    """
    Execute a consolidated antechamber pipeline task inside an isolated sandbox.

    Parameters
    ----------
    coords : list of list of float or numpy.ndarray
        Array (N, 3) tracking the geometry definitions of the system.
    elements : list of str
        Ordered chemical elements tracking index parity exactly.
    charge : int, default 0
        The explicit overall net formal charge of the chemical layout.
    get_charges : bool, default True
        If True, extracts calculated atomic partial charges from output streams.
    get_types : bool, default True
        If True, extracts force field parameter designations.
    get_names : bool, default True
        If True, extracts unique non-clashing structural labels.
    charge_method : str, default "bcc"
        Target calculation matrix approach. Supported options: 'bcc', 'abcg2'.
    atom_type_family : str, default "gaff"
        Target force field parametrization style. Supported:
        "gaff", "amber", "gaff2", "bcc", "abcg2", "sybyl"
    optimize_geometry : bool, default False
        If True, executes semi-empirical geometry relaxation optimizations.
    adjust_names : str, default None
        Optional structural atom label adjustment configuration string (-an).
    base_dir : str, default "tmpfiles"
        The file directory path context targeted for temporary computation data.
    resname : str, default "MOL"
        Optional residue name
    verbose : bool, default False
        If True, diagnostic steps print directly to stdout. If False, output is
        only displayed on execution failure.

    Returns
    -------
    results : dict
        A lookup interface tracking structural metadata arrays. Contains keys:
        "charges": list of float or None
        "types": list of str or None
        "names": list of str or None
    """
    import os
    import io
    import sys
    import subprocess
    import tempfile
    import parmed

    log_buf = io.StringIO()
    exe = _check_antechamber_executable()
    
    if atom_type_family.lower() not in ["gaff", "amber", "gaff2", "bcc", "abcg2", "sybyl"]:
        raise ValueError("Atom type family selection must be 'gaff' or 'amber'.")
        
    clean_method = charge_method.lower()
    if clean_method in ["bcc", "am1bcc"]:
        method_flag = "bcc"
    elif clean_method == "abcg2":
        method_flag = "abcg2"
    else:
        raise ValueError(f"Unsupported charge calculation method: {charge_method}")
        
    try:
        if not os.path.exists(base_dir):
            os.makedirs(base_dir, exist_ok=True)
            
        input_rel = _generate_high_precision_input_mol2(
            coords, elements, base_dir, charge=charge, log_buffer=log_buf, resname=resname
        )
        input_abs = os.path.abspath(input_rel)
        
        out_dir_rel = tempfile.mkdtemp(dir=base_dir)
        out_dir_abs = os.path.abspath(out_dir_rel)
        mol2_out_abs = "/".join([out_dir_abs, "output.mol2"])
        
        cmd = [
            exe,
            "-i", input_abs,
            "-fi", "mol2",
            "-o", mol2_out_abs,
            "-fo", "mol2",
            "-at", atom_type_family.lower(),
            "-seq", "n"
        ]
        
        if get_charges:
            cmd.extend([
                "-c", method_flag,
                "-nc", str(charge)
            ])
            if not optimize_geometry:
                maxcyc = 999 if optimize_geometry else 0
                ek_str = f"qm_theory='AM1', qmcharge={charge}, maxcyc={maxcyc}"
                cmd.extend([
                    "-ek", ek_str
                ])
            
        if adjust_names is not None:
            cmd.extend(["-an", adjust_names])
            
        log_buf.write("--- HYPOTHESIS TESTING DIAGNOSTICS: RUNNING WITH SUBPROCESS CWD ---\n")
        log_buf.write(f"Parent process working directory: {os.getcwd()}\n")
        log_buf.write(f"Child process target cwd workspace sandbox: {out_dir_abs}\n")
        cmdstr = " ".join(cmd)
        log_buf.write(f"Executing compiled command sequence array: {cmdstr}\n")
        
        res = subprocess.run(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True, 
            check=False,
            cwd=out_dir_abs
        )
        
        log_buf.write("--- EXECUTION STDOUT STREAMS ---\n")
        log_buf.write(f"{res.stdout}\n")
        log_buf.write("--- EXECUTION STDERR STREAMS ---\n")
        log_buf.write(f"{res.stderr}\n")
        
        if res.returncode != 0 or not os.path.exists(mol2_out_abs):
            raise RuntimeError("Antechamber consolidated execution failure generated.")
            
        loaded_mol2 = parmed.load_file(mol2_out_abs, structure=True)
        
        if len(loaded_mol2.atoms) != len(elements):
            log_buf.write("DIAGNOSTIC CRITICAL FAIL: Element row alignment breakdown.\n")
            log_buf.write(f"Input nodes size: {len(elements)} vs Generated output: {len(loaded_mol2.atoms)}\n")
            
        raw_charges = [float(atom.charge) for atom in loaded_mol2.atoms] if get_charges else None
        fixed_charges = FixCharges(raw_charges, digits=4) if get_charges else None

        output_dict: Dict[str, Optional[List[Union[str, float]]]] = {
            "charges": fixed_charges,
            "types": [str(atom.type) for atom in loaded_mol2.atoms] if get_types else None,
            "names": [atom.name for atom in loaded_mol2.atoms] if get_names else None
        }
        
        if get_charges:
            log_buf.write(f"DIAGNOSTIC BALANCED NET CHARGE CHECK: Sum = {sum(output_dict['charges']):.4f}\n")
            
        if verbose:
            sys.stdout.write(log_buf.getvalue())
            sys.stdout.flush()
            
        return output_dict

    except Exception as err:
        sys.stderr.write("\n!!! CRITICAL PIPELINE FAILURE INTERCEPTED !!!\n")
        sys.stderr.write("Dumping captured file trace history buffer metrics:\n")
        sys.stderr.write(log_buf.getvalue())
        sys.stderr.write(f"Exception Message Error: {str(err)}\n\n")
        sys.stderr.flush()
        raise


def print_environment_versions():
    """Print system version metrics with stripped terminal color formatting codes."""
    import sys
    import re
    import subprocess
    
    print(f"Python Version: {sys.version}")
    try:
        import parmed
        print(f"ParmEd version: {parmed.__version__}")
    except ImportError:
        print("ParmEd package is not currently installed.")
        
    try:
        from openbabel import pybel
        print("OpenBabel/Pybel chemical perception modules successfully localized.")
    except ImportError:
        print("OpenBabel package is missing from active environment paths.")
        
    try:
        res = subprocess.run(
            ["antechamber", "-h"], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True
        )
        
        ansi_macro = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
        cleaned_help_text = ansi_macro.sub('', res.stdout)
        
        print("\n--- ANTECHAMBER SUBPROCESS HELP SIGNATURE ---")
        for line in cleaned_help_text.split('\n')[:5]:
            if line.strip():
                print(line)
        print("---------------------------------------------\n")
    except Exception as e:
        print(f"Failed to interrogate antechamber directly: {str(e)}")


        


# def run_antechamber_from_structs(
#     structs: List[ffpopt.Struct.Struct],
#     get_charges: bool = True,
#     get_types: bool = True,
#     get_names: bool = True,
#     charge_method: str = "bcc",
#     atom_type_family: str = "gaff",
#     optimize_geometry: bool = False,
#     adjust_names: Optional[str] = None,
#     base_dir: str = "tmpfiles",
#     resname: Optional[str] = "MOL",
#     verbose: bool = False
# ) -> Tuple[ffpopt.Struct.ListOfStruct,ffpopt.Struct.ListOfStruct]:
#     """
#     Execute a consolidated antechamber pipeline task inside an isolated sandbox.

#     Parameters
#     ----------
#     structs : list of ffpopt.Struct.Struct
#         The list of structures/conformations
#     get_charges : bool, default True
#         If True, extracts calculated atomic partial charges from output streams.
#     get_types : bool, default True
#         If True, extracts force field parameter designations.
#     get_names : bool, default True
#         If True, extracts unique non-clashing structural labels.
#     charge_method : str, default "bcc"
#         Target calculation matrix approach. Supported options: 'bcc', 'abcg2'.
#     atom_type_family : str, default "gaff"
#         Target force field parametrization style. Supported:
#         "gaff", "amber", "gaff2", "bcc", "abcg2", "sybyl"
#     optimize_geometry : bool, default False
#         If True, executes semi-empirical geometry relaxation optimizations.
#     adjust_names : str, default None
#         Optional structural atom label adjustment configuration string (-an).
#     base_dir : str, default "tmpfiles"
#         The file directory path context targeted for temporary computation data.
#     resname : str, default "MOL"
#         Optional residue name
#     verbose : bool, default False
#         If True, diagnostic steps print directly to stdout. If False, output is
#         only displayed on execution failure.

#     Returns
#     -------
#     results : tuple of ListOfStruct, ListOfStruct
#         The two ListOfStruct are identical except for the charges.
#         The first ListOfStruct has conformer-specific charges.
#         The second ListOfStruct uses the conformer-averaged charges.
#     """

#     import copy
#     from . Struct import ListOfStruct
    
#     ostructs = copy.deepcopy(structs)
#     for s in ostructs:
#         res = run_antechamber\
#             ( s.data["positions"],
#               s.data["elements"],
#               s.GetCharge(),
#               get_charges=get_charges,
#               get_types=get_types,
#               get_names=get_names,
#               charge_method=charge_method,
#               atom_type_family=atom_type_family,
#               optimize_geometry=optimize_geometry,
#               adjust_names=adjust_names,
#               base_dir=base_dir,
#               resname=resname,
#               verbose=verbose )
#         if get_names:
#             s.data["names"] = res["names"]
#         if get_types:
#             s.data["types"] = res["types"]
#         if get_charges:
#             s.data["charges"] = res["charges"]

#     astructs = copy.deepcopy(ostructs)
#     if get_charges and len(ostructs) > 1:
#         qs = []
#         for s in astructs:
#             qs.append( s.data["charges"] )
#         qs = np.array(qs)
#         qs = qs.mean(qs,dim=0).tolist()
#         qs = FixCharges(qs)
#         for s in astructs:
#             s.data["charges"] = qs

#     return ListOfStruct(ostructs), ListOfStruct(astructs)
        
