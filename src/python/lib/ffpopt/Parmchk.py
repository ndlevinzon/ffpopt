#!/usr/bin/env python3
import subprocess
import os
import sys
import argparse
import parmed as pmd


# Map all standard user string/integer inputs to parmchk2 integer flags and Amber filename keys
FF_MAPPING = {
    "1": ("1", "gaff.dat"), "gaff": ("1", "gaff.dat"),
    "2": ("2", "gaff2.dat"), "gaff2": ("2", "gaff2.dat"),
    "3": ("3", "parm99.dat"), "parm99": ("3", "parm99.dat"),
    "4": ("4", "parm10.dat"), "parm10": ("4", "parm10.dat"),
    "5": ("5", "lipid14.dat"), "lipid14": ("5", "lipid14.dat")
}





def extract_all_mass_and_nonbon(input_mol2, ff_flag):
    """
    Directly replaces extract_all_mass_and_nonbon. Parses the mol2 file and
    the parameter set via ParmEd, completely eliminating parmchk2 and temp files.
    """
    import os
    import parmed as pmd
    from parmed.amber.parameters import AmberParameterSet

    mass_dict = {}
    nonbon_dict = {}

    # 1. Parse the input mol2 file using ParmEd to extract unique atom types
    try:
        mol = pmd.load_file(input_mol2)
    except Exception as e:
        print(f"[ExtractAllMassAndNonbon] Error loading mol2 file: {e}")
        return mass_dict, nonbon_dict

    unique_types = {atom.type for atom in mol.atoms}

    # 2. Extract the parameter file name using the global FF_MAPPING structure
    ff_entry = FF_MAPPING.get(ff_flag) or FF_MAPPING.get(str(ff_flag))
    if not ff_entry:
        print(f"[ExtractAllMassAndNonbon] Error: Force field flag '{ff_flag}' not found in FF_MAPPING.")
        return mass_dict, nonbon_dict
        
    filename = ff_entry[1]  # Extract the file string (e.g., 'parm10.dat') from the tuple
    amber_home = os.environ.get("AMBERHOME", "")
    
    # 3. Locate the parameter file path
    possible_paths = [
        filename,
        os.path.join(amber_home, "dat", "leap", "parm", filename),
        os.path.join(amber_home, "dat", "parm", filename)
    ]
    ff_path = next((path for path in possible_paths if os.path.exists(path)), None)
    if not ff_path:
        print(f"[ExtractAllMassAndNonbon] Error: Could not find parameter file {filename}")
        return mass_dict, nonbon_dict

    # 4. Load the parameter library exactly once
    params = AmberParameterSet(ff_path)

    # 5. Extract masses and non-bonded variables for each atom type present
    for atom_type in unique_types:
        # Resolve Mass (Masses are typically mapped strictly by the literal atom type string)
        if atom_type in params.atom_types:
            at_mass = params.atom_types[atom_type]
            # Formats to standard frcmod MASS block alignment: Type Mass
            mass_dict[atom_type] = f"  {atom_type:<4}  {at_mass.mass:7.3f}"
        else:
            print(f"[ExtractAllMassAndNonbon] Warning: Mass for type '{atom_type}' not found in {filename}")

        # Resolve Non-Bonded Parameters (Handles native equivalencies like NC -> N)
        target_type = atom_type
        if target_type not in params.atom_types:
            #for eq_list in params.nbequiv:
            #    if atom_type in eq_list:
            #        target_type = eq_list[0]  # The first entry serves as the parameter source
            #        break
            for standard_type in params.atom_types:
                if standard_type.lower() == atom_type.lower():
                    target_type = standard_type
                    break

        if target_type in params.atom_types:
            at_nonbon = params.atom_types[target_type]
            # Formats to standard frcmod NONBON block alignment: Type Radius Epsilon
            nonbon_dict[atom_type] = f"  {atom_type:<4}        {at_nonbon.rmin:6.4f}  {at_nonbon.epsilon:6.4f}"
        else:
            print(f"[ExtractAllMassAndNonbon] Warning: Non-bonded parameters for type '{atom_type}' could not be resolved.")

    #print("Read NONBON:", nonbon_dict)
    return mass_dict, nonbon_dict





# Comprehensive, generalized GAFF2/Amber chemical family fallback wildcards
GAFF2_WILDCARDS = {
    # --- CARBONS ---
    'c3': ['c2', 'ca', 'c'], 'cx': ['c3', 'c2', 'c'], 'cy': ['c3', 'c2', 'c'], 'c5': ['c3', 'c2', 'c'], 'c6': ['c3', 'c2', 'c'],
    'c':  ['c2', 'ca', 'c3'], 'cs': ['c', 'c2', 'c3'], 'c2': ['c3', 'ca', 'c'], 'ce': ['c2', 'ca', 'c'], 'cf': ['c2', 'ca', 'c'],
    'cu': ['c2', 'c3', 'c'],  'cv': ['c2', 'c3', 'c'],  'cz': ['c2', 'c3', 'c'],  'ca': ['c2', 'c3', 'c'], 'cp': ['ca', 'c2', 'c'],
    'cq': ['ca', 'c2', 'c'],  'cc': ['ca', 'c2', 'c'],  'cd': ['ca', 'c2', 'c'],  'c1': ['c2', 'c3', 'c'], 'cg': ['c1', 'c2', 'c3'],
    'ch': ['c1', 'c2', 'c3'],
    # --- OXYGENS ---
    'o':  ['o2', 'oh', 'os'], 'o2': ['o', 'oh', 'os'], 'oh': ['os', 'o'], 'os': ['oh', 'o'], 'op': ['oh', 'os', 'o'],
    'oq': ['oh', 'os', 'o'],  'ow': ['oh', 'os'],
    # --- NITROGENS ---
    'n':  ['ns', 'nt', 'n2', 'n3'], 'ns': ['n', 'nt', 'n2', 'n3'], 'nt': ['n', 'ns', 'n2', 'n3'], 'ni': ['n', 'ns', 'n2'],
    'nj': ['n', 'ns', 'n2'],        'n3': ['n4', 'n2', 'n'],        'n4': ['n3', 'n+', 'n'],       'n7': ['n3', 'n2', 'n'],
    'n8': ['n3', 'n2', 'n'],        'n9': ['n3', 'n4', 'n+'],       'n+': ['n4', 'n3', 'n'],       'nx': ['n4', 'n3', 'n'],
    'ny': ['n4', 'n3', 'n'],        'nz': ['n4', 'n3', 'n'],        'nk': ['n4', 'n3', 'n'],       'nl': ['n4', 'n3', 'n'],
    'np': ['n3', 'n2', 'n'],        'nq': ['n3', 'n2', 'n'],        'n5': ['n3', 'n2', 'n'],       'n6': ['n3', 'n2', 'n'],
    'n2': ['n', 'n3', 'na'],        'na': ['n', 'n2', 'nb'],        'nb': ['n', 'n2', 'na'],       'nc': ['n', 'n2', 'nd'],
    'nd': ['n', 'n2', 'nc'],        'ne': ['n', 'n2', 'n3'],       'nf': ['n', 'n2', 'n3'],       'nh': ['n3', 'n2', 'n'],
    'nu': ['nh', 'n3', 'n2'],       'nv': ['nh', 'n3', 'n2'],       'nm': ['nh', 'n3', 'n2'],       'nn': ['nh', 'n3', 'n2'],
    'no': ['n2', 'n', 'n3'],        'n1': ['n2', 'n', 'n3'],
    # --- SULFURS ---
    's':  ['ss', 'sh', 's2', 's4', 's6'], 's2': ['s', 'ss', 'sx'], 's4': ['s', 'ss', 'sx'], 's6': ['s', 'ss', 'sy'],
    'sx': ['s4', 's2', 's'],               'sy': ['s6', 's', 'ss'], 'ss': ['sh', 's'],     'sh': ['ss', 's'],
    'sp': ['ss', 'sh', 's'],               'sq': ['ss', 'sh', 's'],
    # --- PHOSPHORUS ---
    'p':  ['p5', 'p4', 'p3'], 'p3': ['p', 'p5'], 'p2': ['p3', 'p'], 'p4': ['p5', 'px', 'p'], 'p5': ['p4', 'py', 'p'],
    'px': ['p4', 'p5', 'p'],  'py': ['p5', 'p4', 'p'], 'pb': ['p', 'p3'], 'pc': ['p', 'p3'], 'pd': ['p', 'p3'],
    'pe': ['p', 'p3'],        'pf': ['p', 'p3'],
    # --- HALOGENS ---
    'f':  ['cl', 'br', 'i'], 'cl': ['br', 'f', 'i'], 'br': ['cl', 'f', 'i'], 'i':  ['br', 'cl', 'f'],
    # --- HYDROGENS ---
    'hc': ['h1', 'h2', 'h3', 'ha'], 'h1': ['hc', 'h2', 'h3', 'ha'], 'h2': ['hc', 'h3', 'h1', 'ha'], 'h3': ['hc', 'h2', 'h1', 'ha'],
    'h4': ['ha', 'hc', 'h1'],       'h5': ['ha', 'hc', 'h1'],       'ha': ['hc', 'h1', 'hn'],       'hx': ['hc', 'h1', 'ha'],
    'hn': ['h1', 'ha', 'hc'],       'ho': ['h1', 'ha', 'hn', 'hc'], 'hp': ['ho', 'h1', 'hn'],       'hs': ['hc', 'h1', 'ho'],
    'hw': ['ho', 'h1'],             'hb': ['hc', 'h1']
}

def get_mode(line,current_mode):
    line_upper = line.upper()
    if "MASS" in line_upper: current_mode = "MASS"
    elif "BOND" in line_upper: current_mode = "BOND"
    elif "ANGL" in line_upper: current_mode = "ANGLE"
    elif "DIHE" in line_upper: current_mode = "DIHEDRAL"
    elif "IMPR" in line_upper: current_mode = "IMPROPER"
    elif "NONB" in line_upper: current_mode = "NONBON"
    return current_mode

def tokenize(line,current_mode):
    tokens = line.upper().split()
    dash_count = 0
    target_count = 0
    if current_mode == "BOND":
        target_count = 1
    elif current_mode == "ANGLE":
        target_count = 2
    elif current_mode == "DIHEDRAL":
        target_count = 3
    elif current_mode == "IMPROPER":
        target_count = 3

    for i in range(len(tokens)):
        if len(tokens) < 2:
            break
        count = tokens[0].count("-")
        if count >= target_count:
            break
        if count < target_count:
            tokens[0:2] = [ tokens[0] + " " + tokens[1] ]
    return tokens


def get_attn_lines_with_modes(filepath):
    attn_entries = []
    if not os.path.exists(filepath):
        return attn_entries
    current_mode = None
    with open(filepath, 'r') as f:
        for line in f:
            current_mode = get_mode(line,current_mode)
            if "ATTN" in line or "need revision" in line:
                tokens = tokenize(line,current_mode)
                if tokens:
                    attn_entries.append((current_mode, tokens))
    return attn_entries


def find_matching_atom_indices(mode, amber_key, amber_structure):
    target_types = [ a.strip() for a in amber_key.split('-') ]
    rev_types = target_types[::-1]
    if mode == "BOND" and len(target_types) == 2:
        for bond in amber_structure.bonds:
            t1, t2 = bond.atom1.type, bond.atom2.type
            tlist = [t1,t2]
            if tlist == target_types or tlist == rev_types:
                return [bond.atom1.idx + 1, bond.atom2.idx + 1]
    elif mode == "ANGLE" and len(target_types) == 3:
        for b1 in amber_structure.bonds:
            for b2 in amber_structure.bonds:
                if b1 == b2: continue
                b1_atoms = {b1.atom1, b1.atom2}
                b2_atoms = {b2.atom1, b2.atom2}
                shared = b1_atoms.intersection(b2_atoms)
                if len(shared) == 1:
                    #print(shared)
                    center_atom = list(shared)
                    outer1_atom = (b1_atoms - shared).pop()
                    outer2_atom = (b2_atoms - shared).pop()
                    #print(target_types,center_atom[0],outer1_atom,outer2_atom)
                    t1, t2, t3 = outer1_atom.type, center_atom[0].type, outer2_atom.type
                    tlist = [t1,t2,t3]
                    if tlist == target_types or tlist == rev_types:
                        return [outer1_atom.idx + 1, center_atom[0].idx + 1, outer2_atom.idx + 1]
    elif mode == "NONBON" and len(target_types) == 1:
        for a in amber_structure.atoms:
            if a.type == target_types[0]:
                return [ a.idx + 1 ]

    return None
def load_master_parmed_db(ff_filename):
    amber_home = os.environ.get("AMBERHOME")
    if not amber_home:
        print("Error: AMBERHOME is not set.")
        return None
    db_path = os.path.join(amber_home, "dat", "leap", "parm", ff_filename)
    if not os.path.exists(db_path):
        print("Error: Could not find parameter file at:", db_path)
        return None
    return pmd.amber.AmberParameterSet(db_path)


def query_parmed_with_wildcards(ff_params, mode, gaff_types, verbose):
    if mode == "BOND":
        t1, t2 = gaff_types
        for alt_t1 in [t1] + GAFF2_WILDCARDS.get(t1, []):
            for alt_t2 in [t2] + GAFF2_WILDCARDS.get(t2, []):
                key = (min(alt_t1, alt_t2), max(alt_t1, alt_t2))
                if key in ff_params.bond_types:
                    b_type = ff_params.bond_types[key]
                    if verbose:
                        print("[RunParmchk] Found wildcard match:", key)
                    return f"  {b_type.k:8.2f} {b_type.req:8.3f}"
    elif mode == "ANGLE":            
        t1, t2, t3 = gaff_types
        for alt_t2 in [t2] + GAFF2_WILDCARDS.get(t2, []):
            for alt_t1 in [t1] + GAFF2_WILDCARDS.get(t1, []):
                for alt_t3 in [t3] + GAFF2_WILDCARDS.get(t3, []):
                    key = (min(alt_t1, alt_t3), alt_t2, max(alt_t1, alt_t3))
                    
                    if key in ff_params.angle_types:
                        a_type = ff_params.angle_types[key]
                        if verbose:
                            print("[RunParmchk] Found wildcard match:", key)
                        return f"  {a_type.k:8.2f} {a_type.theteq:8.3f}"
    elif mode == "NONBON":
        t1 = gaff_types[0]
        for alt_t1 in [t1] + GAFF2_WILDCARDS.get(t1, []):
            key = t1
            if key in ff_params.atom_types:
                a_type = ff_params.atom_types[key]
                if verbose:
                    print("[RunParmchk] Found wildcard match:", key)
                return f"  {a_type.rmin:8.4f} {a_type.epsilon:8.4f}"

    return None

def run_parmchk2_OLD(input_mol2, output_frcmod, ff_selection="1", ff_frcmod=None, keep_files=False, all_parameters=False, verbose=True ):
    """
    Run the Amber ``parmchk2`` utility to generate a force field modification file.

    This function serves as a Python interface to the ``parmchk2`` command-line
    tool. It extracts mass and non-bonded parameter databases, validates the 
    molecular topology against a chosen force field dataset, and produces a 
    corresponding frcmod file containing missing or modified parameters.

    Parameters
    ----------
    input_mol2 : str
        Path to the input molecular file in Tripos MOL2 format.
    output_frcmod : str
        Path where the output force field modification (frcmod) file will 
        be saved.
    ff_selection : str, default "1"
        The force field parameter set selection identifier. Supported values 
        include:
        
        * ``"1"`` or ``"gaff"`` : General Amber Force Field (GAFF)
        * ``"2"`` or ``"gaff2"`` : General Amber Force Field v2 (GAFF2)
        * ``"3"`` or ``"parm99"`` : Amber parm99
        * ``"4"`` or ``"parm10"`` : Amber parm10
        * ``"5"`` or ``"lipid14"`` : Amber lipid14
    ff_frcmod : str, optional
        A plus-separated string specifying standard Amber frcmod files 
        to be loaded for specific macromolecular contexts. 
        Example: ``"ff14SB+bsc1+yil"``.
    keep_files : bool, default False
        If True, retains intermediate or temporary files generated during 
        the database extraction steps.
    all_parameters : bool, default False
        If True, passes the ``-a Y`` flag to ``parmchk2`` to print out all 
        force field parameters, including those already defined within the 
        base parameter files. If False, prints only the missing parameters.
    verbose : bool, default True
        If True, print debug info to stdout
    
    Returns
    -------
    None

    Raises
    ------
    KeyError
        If the parsed ``ff_selection`` string does not exist in the 
        global ``FF_MAPPING`` dictionary.
    subprocess.CalledProcessError
        If the underlying ``parmchk2`` subprocess execution fails.

    Examples
    --------
    >>> run_parmchk2("ligand.mol2", "ligand.frcmod", ff_selection="gaff2")
    Step 1: Extracting full MASS/NONBON parameters database...
    Step 2: Checking parameters topology context...
    """
    
    ff_flag, ff_file = FF_MAPPING[str(ff_selection).lower()]

    if ff_selection == "gaff":
        ff_selection = "1"
    elif ff_selection == "gaff2":
        ff_selection = "2"
    elif ff_selection == "parm99":
        ff_selection = "3"
    elif ff_selection == "parm10":
        ff_selection = "4"
    elif ff_selection == "lipid14":
        ff_selection = "5"

    if verbose:
        print("[RunParmchk] Step 1: Extracting full MASS/NONBON parameters database...")
    master_mass, master_nonbon = extract_all_mass_and_nonbon(input_mol2, ff_flag)

    if verbose:
        print("[RunParmchk] Step 2: Checking parameters topology context...")
    options = ["parmchk2", "-i", input_mol2, "-f", "mol2", "-o", output_frcmod, "-s", ff_flag]
    if ff_frcmod is not None:
        options.extend( [ "-frc", ff_frcmod ] )
    if all_parameters:
        options.extend( [ "-a", "Y" ] )
    subprocess.run(options)
    
    missing_entries = get_attn_lines_with_modes(output_frcmod)

    if verbose:
        print("[RunParmchk] Step 3: Loading databases...")
    amber_struct = pmd.load_file(input_mol2, structure=True)
    
    # ff_params = load_master_parmed_db(ff_file)
    # if not ff_params:
    #     print("[RunParmchk] Error: Could not load target master force field file:", ff_file)
    #     return False
        
    temp_gaff_mol2 = "temp_gaff2_reference.mol2"
    temp_gaff_frcmod = "temp_gaff2_reference.frcmod"
    
    try:
        fixed_lines = []
        fixed_count = 0
        
        if missing_entries:
            if verbose:
                print("[RunParmchk] Step 4: Creating temporary GAFF2 reference files for missing terms...")
            cmd = ["antechamber", "-i", input_mol2, "-fi", "mol2", "-o", temp_gaff_mol2, "-fo", "mol2", "-at", "gaff2", "-cf", "cc"]
            subprocess.run(cmd, capture_output=True, text=True)
            gaff2_struct = pmd.load_file(temp_gaff_mol2, structure=True)
            subprocess.run(["parmchk2", "-i", temp_gaff_mol2, "-f", "mol2", "-o", temp_gaff_frcmod, "-s", "2", "-a", "Y"], capture_output=True)


            gaff2_ff_flag, gaff2_ff_file = FF_MAPPING["gaff2"]
            ff_params = load_master_parmed_db(gaff2_ff_file)
            if not ff_params:
                print("Error: Could not load target master force field file:", gaff2_ff_file)
                return False

            
            if verbose:
                print("[RunParmchk] Step 5: Patching lines...")
            current_mode = None
            with open(output_frcmod, 'r') as f:
                for line in f:
                    if "REMARK LINE GOES HERE" in line.upper():
                        continue
                    line_upper = line.upper()
                    current_mode = get_mode(line,current_mode)
                    
                    if "ATTN" in line or "need revision" in line:
                        tokens = tokenize(line,current_mode)
                        if tokens:
                            amber_key = tokens[0]
                            if verbose:
                                print("[RunParmchk] Processing line:", line.rstrip())
                            atom_indices = find_matching_atom_indices(current_mode, amber_key, amber_struct)
                            if atom_indices:
                                types_list = [gaff2_struct.atoms[idx-1].type for idx in atom_indices]
                                gaff_types = tuple(types_list)
                                matched_value = query_parmed_with_wildcards(ff_params, current_mode, gaff_types,verbose)
                                
                                if not matched_value:
                                    print("[debug] Failed to find suitable parameters for this line:",line)
                                    print("[debug] tokens",tokens)
                                    print("[debug] atom_indices",atom_indices)
                                    print("[debug] gaff_types",gaff_types)
                                
                                if matched_value:
                                    spacing = " " * (10 - len(str(amber_key)))
                                    line = str(amber_key) + spacing + str(matched_value) + " Fixed via lookup\n"
                                    fixed_count += 1
                    fixed_lines.append(line)
        else:
            if verbose:
                print("[RunParmchk] No missing interactions detected. Proceeding to enforce complete file parameters...")
            with open(output_frcmod, 'r') as f:
                for line in f:
                    if "REMARK LINE GOES HERE" in line.upper():
                        continue
                    fixed_lines.append(line)
                
        seen_atom_types = set()
        if all_parameters:
            seen_atom_types = set(master_mass.keys())
        else:
            prefixes = ("BON", "ANG", "DIH", "IMP", "NON")
            current_mode = None
            for line in fixed_lines:
                line_upper = line.upper()
                current_mode = get_mode(line,current_mode)
                tokens = tokenize(line,current_mode)
                if tokens and current_mode in ["BOND", "ANGLE", "DIHEDRAL", "IMPROPER", "NONBON"]:
                    if tokens[0].startswith(prefixes):
                        continue
                    # Safely split by hyphen on the identifier string itself (tokens[0])
                    for part in tokens[0].split('-'):
                        clean_type = part.strip()
                        if clean_type and not clean_type.replace('.', '', 1).isdigit():
                            #print("add ",clean_type)
                            seen_atom_types.add(clean_type)
        if verbose:
            print(f"[RunParmchk]  Final list of atom types to write: {list(seen_atom_types)}")
        
        final_lines = []
        final_lines.append("Generated by parmchk2 Wrapper Workflow\n")
        
        final_lines.append("MASS\n")
        for at_type in sorted(seen_atom_types):
            if at_type in master_mass:
                final_lines.append(master_mass[at_type] + "\n")
        final_lines.append("\n")
        
        skip_mode = None
        for l in fixed_lines:
            l_upper = l.upper()
            if "MASS" in l_upper: skip_mode = "MASS"; continue
            if "NONB" in l_upper: skip_mode = "NONBON"; continue
            if any(h in l_upper for h in ["BOND", "ANGL", "DIHE", "IMPR"]): skip_mode = None
            if skip_mode == "MASS" or skip_mode == "NONBON": continue
            final_lines.append(l)
            
        final_lines.append("NONBON\n")
        for at_type in sorted(seen_atom_types):
            if at_type in master_nonbon:
                final_lines.append(master_nonbon[at_type] + "\n")
        final_lines.append("\n")
        
        with open(output_frcmod, 'w') as f:
            f.writelines(final_lines)
        if verbose:
            print("[RunParmchk] Step 6: Complete Workflow Modification Done.")
            print("[RunParmchk] Fixed parameters count:", fixed_count)
        
    finally:
        cleanup_targets = ["ANTECHAMBER.ESP", "ANTECHAMBER_AC.AC", "ANTECHAMBER_AC.AC0"]
        if not keep_files:
            cleanup_targets.extend([temp_gaff_mol2, temp_gaff_frcmod])
        for temp_file in cleanup_targets:
            if os.path.exists(temp_file):
                os.remove(temp_file)


#######################################################################################################
#######################################################################################################



def run_parmchk2(input_mol2, output_frcmod, ff_selection="1", ff_frcmod=None, keep_files=False, all_parameters=False, verbose=True, ff_missing="1", score_file=None, ff_afrc=None, only_missing=False ):
    """
    Run the Amber ``parmchk2`` utility to generate a force field modification file.

    This function serves as a Python interface to the ``parmchk2`` command-line
    tool. It extracts mass and non-bonded parameter databases, validates the 
    molecular topology against a chosen force field dataset, and produces a 
    corresponding frcmod file containing missing or modified parameters.

    Parameters
    ----------
    input_mol2 : str
        Path to the input molecular file in Tripos MOL2 format.
    output_frcmod : str
        Path where the output force field modification (frcmod) file will 
        be saved.
    ff_selection : str, default "1"
        The force field parameter set selection identifier. Supported values 
        include:
        
        * ``"1"`` or ``"gaff"`` : General Amber Force Field (GAFF)
        * ``"2"`` or ``"gaff2"`` : General Amber Force Field v2 (GAFF2)
        * ``"3"`` or ``"parm99"`` : Amber parm99
        * ``"4"`` or ``"parm10"`` : Amber parm10
        * ``"5"`` or ``"lipid14"`` : Amber lipid14
    ff_frcmod : str, optional
        A plus-separated string specifying standard Amber frcmod files 
        to be loaded for specific macromolecular contexts. 
        Example: ``"ff14SB+bsc1+yil"``.
    keep_files : bool, default False
        If True, retains intermediate or temporary files generated during 
        the database extraction steps.
    all_parameters : bool, default False
        If True, passes the ``-a Y`` flag to ``parmchk2`` to print out all 
        force field parameters, including those already defined within the 
        base parameter files. If False, prints only the missing parameters.
    verbose : bool, default True
        If True, print debug info to stdout

    ff_missing : str, default "1"
       The force field parameter set identifier used to lookup missing parameters.
       The only supported values are:
    
        * ``"1"`` or ``"gaff"`` : General Amber Force Field (GAFF)
        * ``"2"`` or ``"gaff2"`` : General Amber Force Field v2 (GAFF2)

        The default is to lookup missing parameters from gaff.

    score_file : str, optional (default=None)
       This is the name of a supplemental atom type score file. The file is
       used by parmchk to lookup parameters for similar atom types based on
       a score. The standard set of amber and gaff atom types are in a
       global score file ${AMBERHOME}/dat/antechamber/PARMCHK.DAT.
       If you introduce new atom types, you'll either need to modify that
       file or create additional entries in a separate file. The score_file
       is the separate file. It is directly passed to parmchk2 via the -atc
       option. For more information, see "parmchk2 -l"

    ff_afrc : str, optional (default=None)
       This is the -afrc parameter of parmchk2. It is an additional frcmod
       file that supplements the parameters read by -s and -frc.

    only_missing : bool, default=False
       If True, then the output frcmod file only contains the missing parameters.
       If False (default), it includes MASS and NONBON entries for all types
       appearing within the file.
    
    Returns
    -------
    None

    Raises
    ------
    KeyError
        If the parsed ``ff_selection`` string does not exist in the 
        global ``FF_MAPPING`` dictionary.
    subprocess.CalledProcessError
        If the underlying ``parmchk2`` subprocess execution fails.

    Examples
    --------
    >>> run_parmchk2("ligand.mol2", "ligand.frcmod", ff_selection="gaff2")
    Step 1: Extracting full MASS/NONBON parameters database...
    Step 2: Checking parameters topology context...
    """
    
    from pathlib import Path
    
    ff_flag, ff_file = FF_MAPPING[str(ff_selection).lower()]

    if ff_selection == "gaff":
        ff_selection = "1"
    elif ff_selection == "gaff2":
        ff_selection = "2"
    elif ff_selection == "parm99":
        ff_selection = "3"
    elif ff_selection == "parm10":
        ff_selection = "4"
    elif ff_selection == "lipid14":
        ff_selection = "5"

    options = ["parmchk2", "-i", input_mol2, "-f", "mol2", "-o", output_frcmod, "-s", ff_flag]
    if ff_frcmod is not None:
        options.extend( [ "-frc", ff_frcmod ] )
    if all_parameters:
        options.extend( [ "-a", "Y" ] )

    if score_file is not None:
        p = Path(score_file)
        if not p.is_file:
            raise Exception(f"score_file does not exist: {score_file}")
        options.extend( [ "-atc", score_file ] )

    if ff_afrc is not None:
        if not Path(ff_afrc).is_file:
            print(f"[RunParmchk] -afrc {ff_afrc} file not found")
        else:
            options.extend( [ "-afrc", ff_afrc ] )
        
    if verbose:
        print("[RunParmchk] Running: %s"%(" ".join(options)))
        
    subprocess.run(options)
    
    missing_entries = get_attn_lines_with_modes(output_frcmod)

    if verbose:
        print("[RunParmchk] Looking for missing parameters...")
        if len(missing_entries) == 0:
            print("[RunParmchk] No missing interactions detected. Proceeding to enforce complete file parameters...")
        else:
            print(f"[RunParmchk] Number of lines with 'ATTN': {len(missing_entries)}")

            
    out_frcmod = Path(output_frcmod)
    if len(missing_entries) > 0:

        if ff_missing == "gaff":
            ff_missing = "1"
        elif ff_missing == "gaff2":
            ff_missing = "2"
        if not (ff_missing == "1" or ff_missing == "2"):
            raise Exception(f"Invalid ff_missing parameter {ff_missing}. Expected either '1', '2', 'gaff', or 'gaff2'")
        
        tmp_frcmod = out_frcmod.with_suffix(out_frcmod.suffix + ".tmp")
        options = ["parmchk2",
                   "-i", str(out_frcmod),
                   "-f", "frcmod",
                   "-o", str(tmp_frcmod),
                   "-s", ff_missing,
                   "-att", "2"]
        if score_file is not None:
            options.extend( [ "-atc", score_file ] )


        if verbose:
            print("[RunParmchk] Running: %s"%(" ".join(options)))
        subprocess.run(options)

        if tmp_frcmod.exists():
            tmp_frcmod.rename(out_frcmod)
            if verbose:
                print(f"[RunParmchk] Renamed: mv {tmp_frcmod} {out_frcmod}")
        else:
            print(f"[FunParmchk] Failed to create {tmp_frcmod}")

        missing_entries = get_attn_lines_with_modes(output_frcmod)
        if verbose:
            print(f"[RunParmchk] Rechecking the number of lines with 'ATTN': {len(missing_entries)}")
        if len(missing_entries) > 0:
            print(f"[RunParmchk] WARNING: There are still {len(missing_entries)} lines with 'ATTN'")



            
    #fixed_lines = out_frcmod.read_text().splitlines()
    fixed_lines = []
    with open(output_frcmod, 'r') as f:
        for line in f:
            if "REMARK LINE GOES HERE" in line.upper():
                continue
            fixed_lines.append(line)

    mass_lines = {}
    nonb_lines = {}
            
    try:
        seen_atom_types = set()
        if all_parameters:
            seen_atom_types = set(master_mass.keys())
        else:
            prefixes = ("BON", "ANG", "DIH", "IMP", "NON","MAS")
            current_mode = None
            for line in fixed_lines:
                line_upper = line.upper()
                current_mode = get_mode(line,current_mode)
                tokens = tokenize(line,current_mode)
                if tokens and current_mode in "MASS":
                    mass_lines[tokens[0]] = line
                if tokens and current_mode in "NONBON":
                    nonb_lines[tokens[0]] = line
                if tokens and current_mode in ["BOND", "ANGLE", "DIHEDRAL", "IMPROPER", "NONBON"]:
                    if tokens[0].startswith(prefixes):
                        continue
                    # Safely split by hyphen on the identifier string itself (tokens[0])
                    for part in tokens[0].split('-'):
                        clean_type = part.strip()
                        if clean_type and not clean_type.replace('.', '', 1).isdigit():
                            seen_atom_types.add(clean_type)
        if verbose:
            print(f"[RunParmchk] Final list of atom types to write: {list(seen_atom_types)}")

        missing_atom_types = []
        for m in seen_atom_types:
            if m not in mass_lines or m not in nonb_lines:
                missing_atom_types.append( m )

        if verbose:
            print(f"[RunParmchk] Number of missing MASS or NONB types: {len(missing_atom_types)}")


        if not only_missing:
            if len(missing_atom_types) > 0:
                if verbose:
                    print("[RunParmchk] Extracting full MASS/NONBON parameters database...")
                master_mass, master_nonbon = extract_all_mass_and_nonbon(input_mol2, ff_flag)
                for m in seen_atom_types:
                    if m not in mass_lines:
                        if m in master_mass:
                            mass_lines[m] = master_mass[m] + "\n"
                    if m not in nonb_lines:
                        if m in master_nonbon:
                            nonb_lines[m] = master_nonbon[m] + "\n"
            


            
        final_lines = []
        final_lines.append("Generated by parmchk2 Wrapper Workflow\n")
        
        final_lines.append("MASS\n")
        for at_type in sorted(seen_atom_types):
            #if at_type in master_mass:
            #    final_lines.append(master_mass[at_type] + "\n")
            if at_type in mass_lines:
                final_lines.append(mass_lines[at_type])
            else:
                print(f"[RunParmchk] Failed to write MASS entry for {at_type}")
        final_lines.append("\n")
        
        skip_mode = None
        for l in fixed_lines:
            l_upper = l.upper()
            if "MASS" in l_upper: skip_mode = "MASS"; continue
            if "NONB" in l_upper: skip_mode = "NONBON"; continue
            if any(h in l_upper for h in ["BOND", "ANGL", "DIHE", "IMPR"]): skip_mode = None
            if skip_mode == "MASS" or skip_mode == "NONBON": continue
            final_lines.append(l)
            
        final_lines.append("NONBON\n")
        for at_type in sorted(seen_atom_types):
            #if at_type in master_nonbon:
            #    final_lines.append(master_nonbon[at_type] + "\n")
            if at_type in nonb_lines:
                final_lines.append( nonb_lines[at_type] )
            else:
                print(f"[RunParmchk] Failed to write NONB entry for {at_type}")
        final_lines.append("\n")
        
        with open(output_frcmod, 'w') as f:
            f.writelines(final_lines)
        if verbose:
            print("[RunParmchk] Workflow Modification Completed")
        
    finally:
        cleanup_targets = ["ANTECHAMBER.ESP", "ANTECHAMBER_AC.AC", "ANTECHAMBER_AC.AC0"]
        #if not keep_files:
        #    cleanup_targets.extend([temp_gaff_mol2, temp_gaff_frcmod])
        for temp_file in cleanup_targets:
            if os.path.exists(temp_file):
                os.remove(temp_file)







                
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="General parmchk2 wrapper script.")
    parser.add_argument('input_mol2', help="Input mol2 path.")
    parser.add_argument('output_frcmod', help="Output frcmod path.")
    parser.add_argument('-s', default='1', help="Force field parameter selection option flag.")
    parser.add_argument('--keep', action='store_true', help="Keep reference validation files.")
    parser.add_argument('--all', action='store_true', help="Always output all molecule parameters.")
    parser.add_argument('--frc', type=str, required=False, help="frcmod files to be loaded, the supported frcmods include ff99SB, ff14SB, ff03 for proteins , bsc1, ol15, ol3 for DNA and yil for RNA eg. ff14SB+bsc1+yil, ff99SB+bsc1")
    args = parser.parse_args()
    
    if not os.path.exists(args.input_mol2) or str(args.s).lower() not in FF_MAPPING:
        sys.exit(1)
        
    run_parmchk2(args.input_mol2, args.output_frcmod, ff_selection=args.s, ff_frcmod=args.frc, keep_files=args.keep, all_parameters=args.all)
