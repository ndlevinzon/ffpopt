from contextlib import contextmanager, redirect_stderr

@contextmanager
def silence_stderr():
    """Completely silences Python and low-level system subprocess stderr."""
    import os
    # 1. Open a pointer to the OS null device
    null_fd = os.open(os.devnull, os.O_WRONLY)
    # 2. Keep a backup of the real stderr stream
    old_stderr_fd = os.dup(2)
    try:
        # 3. Mute Python's stderr stream layer
        with open(os.devnull, 'w') as f_null:
            with redirect_stderr(f_null):
                # 4. Force lower C/Java system file descriptor 2 to /dev/null
                os.dup2(null_fd, 2)
                yield
    finally:
        # 5. Clean up and restore original stream when exiting the block
        os.dup2(old_stderr_fd, 2)
        os.close(null_fd)
        os.close(old_stderr_fd)



def parse_filename_charge(filename):
    """
    Check if a filename ends with a colon followed by an integer charge, 
    and split it into a text prefix and an integer charge value.

    This function parses file strings that append electrostatic charge metadata 
    directly to the file path (e.g., 'protein.pdb:-2' or 'ligand.mol2:+1'). It 
    safely evaluates leading positive or negative signs on the trailing 
    integer and isolates the core file path prefix.
    
    If the filename does not follow this format (no colon or no trailing 
    integer), the function defaults to treating the entire input string as 
    the prefix with a baseline charge of 0.

    Parameters
    ----------
    filename : str
        The full string or file path to evaluate. Expected format is 
        'path/to/file:charge', where charge is a valid positive or negative 
        integer.

    Returns
    -------
    prefix : str
        The base filename string extracted before the final colon if a match 
        is found; otherwise, the original unaltered `filename` input string.
    charge : int
        The signed integer parsed from the trailing characters if a match 
        is found; otherwise, a default fallback integer value of `0`.

    See Also
    --------
    re.match : Core regular expression matching function utilized.

    Examples
    --------
    >>> parse_filename_charge("molecule.pdb:-1")
    ('molecule.pdb', -1)

    >>> parse_filename_charge("ligand.mol2:+2")
    ('ligand.mol2', 2)

    >>> parse_filename_charge("protein.pdb:3")
    ('protein.pdb', 3)

    >>> parse_filename_charge("my_file_without_charge.pdb")
    ('my_file_without_charge.pdb', 0)

    >>> parse_filename_charge("invalid_format:abc")
    ('invalid_format:abc', 0)
    """
    # Encapsulated import inside the function scope to isolate dependencies
    import re

    # Pattern breakdown:
    # ^(.+)       - Group 1 (prefix): Non-greedily captures all leading characters
    # :           - Matches the literal colon separating the prefix and charge
    # ([+-]?\d+)  - Group 2 (charge): Captures an optional sign (+/-) and digits
    # $           - Assures the pattern matches completely up to the string end
    pattern = r"^(.+):([+-]?\d+)$"
    
    match = re.match(pattern, filename)
    
    if match:
        prefix = match.group(1)
        # Python's int() natively converts signed string integers like "+2" or "-1"
        charge = int(match.group(2))
        return prefix, charge
        
    # Default fallback: Treat the input as the prefix and assign a neutral 0 charge
    return filename, 0



def encode_inchi_urlsafe(inchi_string: str) -> str:
    """Compresses an InChI string and encodes it to URL-safe Base64 (without padding)."""
    import base64
    #import zlib
    # 1. Convert text to raw UTF-8 bytes
    #raw_bytes = inchi_string.encode('utf-8')
    # 2. Compress the bytes losslessly at maximum level
    #compressed_bytes = zlib.compress(raw_bytes, level=9)
    # 3. Encode to URL-safe Base64 (substitutes + and / with - and _)
    #b64_bytes = base64.urlsafe_b64encode(compressed_bytes)
    #b64_string = b64_bytes.decode('ascii')
    # 4. Strip trailing '=' padding characters for cleaner URLs
    #return b64_string.rstrip('=')
    return base64.urlsafe_b64encode(inchi_string.encode()).decode().rstrip('=')


def decode_inchi_urlsafe(encoded_string: str) -> str:
    """Reverses the URL-safe process to completely restore the original InChI string."""
    import base64
    #import zlib
    # 1. Add back the stripped padding ('=') based on string length remainder
    #missing_padding = len(encoded_string) % 4
    #if missing_padding:
    #    encoded_string += '=' * (4 - missing_padding)
    # 2. Decode the URL-safe Base64 string back to compressed binary bytes
    #compressed_bytes = base64.urlsafe_b64decode(encoded_string.encode('ascii'))
    # 3. Decompress the binary bytes back to raw UTF-8 bytes
    #raw_bytes = zlib.decompress(compressed_bytes)
    # 4. Decode bytes back to the original text string
    #return raw_bytes.decode('utf-8')
    encoded_string = encoded_string.rstrip("/").split("/")[-1]
    return base64.urlsafe_b64decode(encoded_string + '=' * (-len(encoded_string) % 4)).decode()






def __GetRDKitMolFromString_backend(fnameormol, quiet=True):
    from rdkit.Chem.inchi import MolFromInchi
    from rdkit.Chem import SanitizeFlags
    from rdkit.Chem import SanitizeMol
    from rdkit.Chem import AddHs
    from rdkit.Chem import MolFromSmiles
    from rdkit import rdBase, Chem

    rdBase.DisableLog('rdApp.warning')
    rdBase.DisableLog('rdApp.error')

    mol = None
    try:
        if not quiet:
            print(f"Trying to interpret {fnameormol} as an inchi string...")
        mol = MolFromInchi(
            fnameormol, removeHs=False, sanitize=False,
            treatWarningAsError=True
        )
        # Sanitize baseline features before adding hydrogens to fix ring networks
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
        mol = AddHs(mol)
        if not quiet:
            print("Success!")
    except:
        if not quiet:
            print("...Failed")
            print(f"Trying to interpret {fnameormol} as a smiles string...")
            
        # CRITICAL FIX: Allow RDKit to sanitize the initial SMILES input string natively.
        # This handles ring counts and double bonds properly before AddHs freezes the topology.
        mol = Chem.MolFromSmiles(fnameormol)
        if mol is not None:
            mol = AddHs(mol)
            if not quiet:
                print("Success!")
        else:
            # Fallback template parse if string rules are highly non-standard
            params = Chem.SmilesParserParams()
            params.sanitize = False
            cansmi = Chem.CanonSmiles(fnameormol)
            mol = Chem.MolFromSmiles(cansmi, params)
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol)
            mol = AddHs(mol)
        
    return mol



def __GetRDKitMolFromString_backend_OLD(fnameormol,quiet=True):
    from rdkit.Chem.inchi import MolFromInchi
    from rdkit.Chem import SanitizeFlags
    from rdkit.Chem import SanitizeMol
    from rdkit.Chem import AddHs
    from rdkit.Chem import MolFromSmiles
    from rdkit.Chem import MolFromMol2Block
    from rdkit import rdBase,Chem
    import parmed
    from io import StringIO

    rdBase.DisableLog('rdApp.warning')
    rdBase.DisableLog('rdApp.error')

    mol = None
    #################################
    try:
        if not quiet:
            print(f"Trying to interpret {fnameormol} as an inchi string...")
        mol = MolFromInchi\
            (fnameormol, removeHs=False, sanitize=False,
             treatWarningAsError=True)
        mol = AddHs(mol)
        #Chem.SanitizeMol(mol)
        if not quiet:
            print("Success!")
    except:
        if not quiet:
            print("...Failed")
            print(f"Trying to interpret {fnameormol} as a smiles string...")
        params = Chem.SmilesParserParams()
        params.sanitize = False
            
        cansmi = Chem.CanonSmiles(fnameormol)
        mol = Chem.MolFromSmiles(cansmi, params)
        mol = AddHs(mol)
        #Chem.SanitizeMol(mol)
        if not quiet:
            print("Success!")
    
    return mol


# def __GetRDKitMolFromString_backend(fnameormol, quiet=False): # Changed quiet=False to see prints
#     from rdkit.Chem.inchi import MolFromInchi
#     from rdkit.Chem import SanitizeFlags, SanitizeMol, AddHs, MolFromSmiles, CanonSmiles
#     from rdkit import rdBase, Chem

#     #print("\n=== STARTING PARSE DEBUGGING ===")
#     #print(f"Input string: {fnameormol}")

#     mol = None
#     try:
#         mol = MolFromInchi(fnameormol, removeHs=False, sanitize=False, treatWarningAsError=True)
#         #print(f"[InChI Step 1] Raw Parse Charge: {Chem.GetFormalCharge(mol)}")
#         mol = AddHs(mol)
#         #print(f"[InChI Step 2] After AddHs Charge: {Chem.GetFormalCharge(mol)}")
#         SanitizeMol(mol)
#         #print(f"[InChI Step 3] After Sanitize Charge: {Chem.GetFormalCharge(mol)}")
#     except Exception as e:
#         #print(f"-> InChI pathway failed/skipped due to: {e}")
#         #print("-> Switching to SMILES pathway...")
        
#         params = Chem.SmilesParserParams()
#         params.sanitize = False
            
#         cansmi = CanonSmiles(fnameormol)
#         mol = MolFromSmiles(cansmi, params)
#         #print(f"[SMILES Step 1] Raw Parse Charge: {Chem.GetFormalCharge(mol)}")
        
#         mol = AddHs(mol)
#         #print(f"[SMILES Step 2] After AddHs Charge: {Chem.GetFormalCharge(mol)}")
        
#         SanitizeMol(mol)
#         #print(f"[SMILES Step 3] After Sanitize Charge: {Chem.GetFormalCharge(mol)}")
    
#     #print(f"=== FINAL RETURNED MOL CHARGE: {Chem.GetFormalCharge(mol)} ===\n")
#     return mol




    
def GetRDKitMolFromString(fnameormol,quiet=True):
    #return __GetRDKitMolFromString_backend(fnameormol,quiet=quiet)
    mol = None
    try:
        mol = __GetRDKitMolFromString_backend(fnameormol,quiet=quiet)
    except:
        fnameormol = decode_inchi_urlsafe(fnameormol)
        mol = __GetRDKitMolFromString_backend(fnameormol,quiet=quiet)
    return mol
    


def RDKitMol2Inchi(mol,tautdep=True):
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import inchi
    options=""
    if tautdep:
        options="/FixedH"
    return inchi.MolToInchi(mol, options=options)


def RDKitMol2InchiToken(mol,tautdep=True):
    s = RDKitMol2Inchi(mol,tautdep=tautdep)
    return encode_inchi_urlsafe(s)
    

def RDKitMol2InchiKey(mol,tautdep=True):
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import inchi
    s = RDKitMol2Inchi(mol,tautdep=tautdep)
    return inchi.InchiToInchiKey(s)


        
# def RDKitMol2Smiles(mol,tautdep=True):
#     import rdkit
#     from rdkit import Chem
#     from rdkit.Chem.MolStandardize import rdMolStandardize
#     if not tautdep:
#         enumerator = rdMolStandardize.TautomerEnumerator()
#         clean_mol = Chem.RemoveHs(mol)
#         tmol = enumerator.Canonicalize(clean_mol)
#     else:
#         tmol = mol
#     return Chem.MolToSmiles(tmol)

def RDKitMol2Smiles(mol, tautdep=False):
    """
    Standalone helper function to convert an RDKit Mol object into a clean SMILES string.
    """
    from rdkit import Chem

    # Make a copy to avoid altering your core GetRDKitAtoms object downstream
    working_mol = Chem.Mol(mol)

    if not tautdep:
        # 1. First, explicitly compute and lock the chiral configurations 
        # while all the explicit hydrogens are physically present in the graph.
        Chem.AssignStereochemistry(working_mol, cleanIt=True, force=True)

        # 2. FIXED: Use the modern RDKit parameters wrapper to strip hydrogens safely.
        # This native object defaults to updating chiral centers automatically.
        params = Chem.RemoveHsParameters()
        working_mol = Chem.RemoveHs(working_mol, params)

        # 3. Re-verify the flags on the newly generated heavy-atom framework
        Chem.AssignStereochemistry(working_mol, cleanIt=True, force=True)

    # 4. Generate the clean, standard canonical SMILES string
    return Chem.MolToSmiles(working_mol)


        


def write_amber_restart(coords, fh, title="Generated by Python"):
    """
    Writes a (Natom, 3) numpy array to an Amber ASCII restart file.
    Format: 2 lines of header followed by coordinates in 6F12.7 format.
    """
    n_atoms = coords.shape[0]
    
    # Flatten coordinates to a 1D array (X1, Y1, Z1, X2, Y2, Z2, ...)
    flat_coords = coords.flatten()

    if True:
        # Line 1: Title (limit to 80 chars)
        fh.write(f"{title[:80]}\n")
        
        # Line 2: Number of atoms (right-justified) and optional time
        # Standard format often uses I5 or I6 for atom count
        fh.write(f"{n_atoms:5d}\n")
        
        # Lines 3+: Coordinates in 6F12.7 format (6 numbers per line, 12 chars each)
        for i in range(0, len(flat_coords), 6):
            chunk = flat_coords[i:i+6]
            line = "".join(f"{val:12.7f}" for val in chunk)
            fh.write(line + "\n")


            

class Struct(object):
    """
    Unified molecular data structure container and interface.

    This class serves as a central data hub that encapsulates structural, 
    topological, and energetic details of a molecular system. It provides 
    seamless factory methods to import structures from diverse backends (such 
    as ParmEd, RDKit, ASE, or file systems) and exposes processing routines 
    tailored for molecular modeling pipelines.

    Attributes
    ----------
    data : dict
        The underlying dictionary housing molecular parameters including 
        atomic positions, element types, partial charges, and bonding graphs.
    restraints : RestraintList or None
        An active list of structural restraints applied to the molecule.
    constraints : ConstraintList or None
        An active list of geometric constraints applied to the system.
    graph : ffpopt.GraphSearch or None
        A search graph representing the covalent connectivity network.
    """
    
    def __init__(self,data):
        """
        Initialize a Struct instance from a raw data configuration dictionary.

        Parameters
        ----------
        data : dict
            A dictionary containing core structural properties. Expected keys 
            include 'positions', 'elements', 'charges', and optionally 
            'resnames', 'residxs', 'constraints', and 'restraints'.
        """
        from copy import deepcopy
        from . Restraints import RestraintList
        from . Constraints import ConstraintList
        
        self.data = deepcopy( data )
        self.restraints = None
        if "restraints" in self.data:
            if self.data["restraints"] is not None:
                self.restraints = RestraintList.from_list_of_dict\
                    ( self.data["restraints"] )
        self.constraints = None
        if "constraints" in self.data:
            if self.data["constraints"] is not None:
                self.constraints = ConstraintList.from_list_of_dict\
                    ( self.data["constraints"] )
                
        self.graph = None

        if "resnames" in self.data:
            resnames = {int(k): v for k, v in self.data["resnames"].items()}
            self.data["resnames"] = resnames
        else:
            self.data["resnames"] = {0: "MOL"}
            
        if "residxs" not in self.data:
            self.data["residxs"] = [0]*len(self.data["elements"])


    @classmethod
    def from_xyz(cls,fname,charge):
        """
        Instantiate a Struct object from a xyz file and charge.

        Parameters
        ----------
        fname : str
            XYZ filename
        charge : int
            Integer net charge

        Returns
        -------
        Struct
            A newly initialized Struct object holding translated properties.
        """
        from . Reader import ReadXyz
        import ase.io
        from collections import defaultdict as ddict
        from . constants   import GetAtomicSymbol, GetAtomicNumber, GetAtomicMass

        netcharge = charge
        mult = 1
        tmpparm = ReadXyz(fname,charge)
        
        bonds = []
        for x in tmpparm.bonds:
            bonds.append( [x.atom1.idx,x.atom2.idx] )

        atypes = []
        for a in tmpparm.atoms:
            atypes.append( a.type )

        atoms = ase.io.read(fname,index=0)
        if "charge" in atoms.info:
            netcharge = atoms.info["charge"]
        if "spin" in atoms.info:
            mult = atoms.info["spin"]
        eles = atoms.get_chemical_symbols()
        crds = atoms.get_positions().tolist()

        netcharge = charge
        
        qs = atoms.get_initial_charges().tolist()
        qsum = sum(qs)
        if qsum != netcharge:
            dq = (netcharge-qsum)/len(eles)
            qs = [ q+dq for q in qs ]

        residxs = [0]*len(eles)
        resnames = {0:"MOL"}
        
        if qs is None:
            qs = [ netcharge / len(eles) ]*len(eles)

        seen = ddict(int)
        names = []
        for e in eles:
            seen[e] += 1
            if seen[e] > 1:
                names.append("%s%i"%(e,seen[e]))
            else:
                names.append(e)

        seen = ddict(int)
        for iname in range(len(names)):
            name = names[iname]
            if name in seen:
                for i in range(98):
                    j=i+1
                    tname = "%s%i"%(name,j)
                    if tname in seen or tname in names:
                        continue
                    else:
                        break
                seen[tname] += 1
                names[iname] = tname

        zs = []
        ms = []
        for e in eles:
            z = GetAtomicNumber(e)
            zs.append(z)
            m = GetAtomicMass(z)
            ms.append(m)

        data = {}
        data["elements"] = eles
        data["types"] = atypes
        data["names"] = names
        data["residxs"] = residxs
        data["resnames"] = resnames
        data["positions"] = crds
        data["charges"] = qs
        data["spin"]  = int(mult)
        data["bonds"] = bonds
        data["parm"] = None
        data["constraints"] = []
        data["restraints"] = []
        data["name"] = "s000"
        data["energy"] = None
        data["forces"] = None
        
        return cls(data)

    

    @classmethod
    def from_parmed(cls,mol):
        """
        Instantiate a Struct object from a ParmEd Structure container.

        Parameters
        ----------
        mol : parmed.Structure
            An active ParmEd molecule instance containing atoms and topologies.

        Returns
        -------
        Struct
            A newly initialized Struct object holding translated properties.
        """
        from . constants import GetAtomicSymbol
        from collections import defaultdict as ddict
        data = {}
        data["positions"] = []
        data["names"] = []
        data["types"] = []
        data["elements"] = []
        data["charges"] = []
        data["residxs"] = []
        resnames = ddict(str)
        q = 0
        for a in mol.atoms:
            data["positions"].append( [a.xx,a.xy,a.xz] )
            data["names"].append(a.name)
            data["types"].append(a.type)
            data["elements"].append( GetAtomicSymbol(a.element) )
            data["charges"].append(a.charge)
            data["residxs"].append(a.residue.idx)
            resnames[a.residue.idx] = a.residue.name
            q += a.charge
        data["resnames"] = dict(resnames)
        data["spin"] = 1
        #data["charge"] = int(round(q))
        data["bonds"] = []
        for x in mol.bonds:
            data["bonds"].append( [x.atom1.idx, x.atom2.idx ] )
        data["parm"] = None
        data["constraints"] = []
        data["restraints"] = []
        data["name"] = "s000"
        data["energy"] = None
        data["forces"] = None
        
        from collections import defaultdict as ddict
        names = data["names"]
        seen = ddict(int)
        for iname in range(len(names)):
            name = names[iname]
            if name in seen:
                for i in range(98):
                    j=i+1
                    tname = "%s%i"%(name,j)
                    if tname in seen or tname in names:
                        continue
                    else:
                        break
                seen[tname] += 1
                names[iname] = tname
        data["names"] = names
        
        return cls(data)

    
    @classmethod
    def from_mol2(cls,fname):
        """
        Instantiate a Struct object by reading a Tripos MOL2 file from disk.

        Parameters
        ----------
        fname : str
            The absolute or relative file path targeting the input MOL2 file.

        Returns
        -------
        Struct
            A newly initialized Struct object containing parsed parameters.
        """
        from . Reader import ReadMol2
        import parmed

        f,q = parse_filename_charge(fname)
        
        mol = ReadMol2(f)

        if f.lower().endswith(".pdb"):
            dq = q / len(mol.atoms)
            for a in mol.atoms:
                a.charge = dq
        
        for a in mol.atoms:
            useelem = isinstance(a.type,int)
            if not useelem:
                useelem = len(a.type) == 0
            if useelem:
                for elem,num in parmed.periodic_table.AtomicNum.items():
                    if num == a.atomic_number:
                        a.type = elem
                        break
        return cls.from_parmed(mol)

    
    @classmethod
    def from_parm(cls,parm7,rst7):
        """
        Instantiate a Struct object from an AMBER topology and coordinate pair.

        Parameters
        ----------
        parm7 : str
            The file path targeting the AMBER parm7/prmtop topology file.
        rst7 : str
            The file path targeting the AMBER rst7/inpcrd coordinate file.

        Returns
        -------
        Struct
            A newly initialized Struct object mapping the parsed AMBER structure.
        """
        import parmed
        mol = parmed.load_file(parm7,xyz=rst7)
        return cls.from_parmed(mol)

    
    @classmethod
    def from_rdkit(cls,mol):
        """
        Instantiate a Struct object from an active RDKit Mol topology object.

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            An active RDKit molecule object container.

        Returns
        -------
        Struct
            A newly initialized Struct object containing translated metrics.
        """
        import numpy as np
        from collections import defaultdict as ddict
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        eles = [atom.GetSymbol() for atom in mol.GetAtoms()]
        try:
            conf = mol.GetConformer()
        except:
            status = AllChem.EmbedMolecule(mol, randomSeed=42)
            conf = mol.GetConformer()

        coords = conf.GetPositions()
        net_charge = Chem.GetFormalCharge(mol)
        #print("[from_rdkit] net_charge: ",net_charge)
        qs = [net_charge/len(eles)]*len(eles)
        bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
        data = {}
        data["positions"] = coords.tolist()
        data["names"] = eles
        data["types"] = eles
        data["elements"] = eles
        data["charges"] = qs

        resnames = ddict(str)
        residxs = []
        for atom in mol.GetAtoms():
            res_info = atom.GetPDBResidueInfo()
            if res_info:
                idx = res_info.GetResidueNumber() - 1
                residxs.append( idx )
                resnames[idx] = res_info.GetResidueName()
            else:
                residxs.append(0)
                resnames[0] = "MOL"
        data["residxs"] = residxs
        data["resnames"] = dict(resnames)
        data["spin"] = 1
        #data["charge"] = int(round(net_charge))
        data["bonds"] = bonds
        data["parm"] = None
        data["constraints"] = []
        data["restraints"] = []
        data["name"] = "s000"
        data["energy"] = None
        data["forces"] = None

        return cls(data)


    def GetCrds(self):
        """
        Return a copy of the atomic coordinates matrix.

        Returns
        -------
        coords : numpy.ndarray
            An (N, 3) matrix mapping precise Cartesian coordinates in Angstroms.
        """
        import numpy as np
        return np.array(self.data["positions"])

    
    def GetAtomicNumbers(self):
        """
        Return the integer atomic numbers of all atoms in the system.

        Returns
        -------
        atomic_numbers : list of int
            An ordered sequence of atomic identifiers matching elements rows.
        """
        from . constants import GetAtomicNumber
        return [ GetAtomicNumber(e) for e in self.data["elements"] ]

    
    def GetCharge(self):
        """
        Compute the absolute rounded net formal charge of the molecule.

        Returns
        -------
        net_charge : int
            The total formal integer net sum of atomic partial charges.
        """
        return int(round(sum(self.data["charges"])))

    
    def ReadAmberParm(self):
        """
        Read an AMBER parm7 topology file from disk and synchronize positions.

        Returns
        -------
        mol : parmed.Structure
            An active ParmEd Structure synchronized with internal positions.

        Raises
        ------
        Exception
            If the 'parm' file tracking key is missing, None, or nonexistent.
        """
        
        from pathlib import Path
        import parmed
        import os
        from tempfile import mkstemp

        name = self.data["name"]
        if "parm" not in self.data:
            raise Exception(f"Key 'parm' not in structure '{name}'")
        pfile = self.data["parm"]
        if pfile is None:
            raise Exception(f"Key 'parm' is None in structure '{name}'")
        if not Path(pfile).exists:
            raise Exception(f"parm7 file {pfile} in structure '{name}' does not exist")
        
        p = parmed.load_file(pfile,structure=True)
        crds = self.data["positions"]
        for i,a in enumerate(p.atoms):
            a.xx = crds[i][0]
            a.xy = crds[i][1]
            a.xz = crds[i][2]
        return p
        

    def GetParmedAtoms(self):
        """
        Generate a fully bonded ParmEd Structure instance of the molecule.

        Returns
        -------
        mol : parmed.Structure
            A newly initialized ParmEd structure container reflecting attributes.
        """
        import parmed as pmd
        from parmed.structure import Structure
        from parmed.topologyobjects import Atom, Bond, Residue
        from . constants import GetAtomicNumber

        
        eles  = self.data["elements"]
        types = self.data["types"]
        qs    = self.data["charges"]
        #q     = self.data["charge"]
        q     = self.GetCharge()
        crds  = self.data["positions"]
        bonds = self.data["bonds"]
        
        # Pull your residue arrays
        if "residxs" in self.data:
            res_indices = self.data["residxs"]    # Array of integer residue indices per atom
        else:
            res_indices = [0] * len(eles)
            
        if "resnames" in self.data:
            res_names_dict = self.data["resnames"]  # Dict mapping int indices to str names
        else:
            res_names_dict = { 0: "MOL" }
        
        if "names" in self.data:
            names = self.data["names"]
        else:
            names = eles

        mol = Structure()
        alist = []
        for i in range(len(eles)):
            z = GetAtomicNumber(eles[i])
            a = Atom(name=names[i],type=types[i],
                     atomic_number=z,charge=qs[i])
            a.xx = crds[i][0]
            a.xy = crds[i][1]
            a.xz = crds[i][2]
            alist.append(a)

            # 1. Grab the specific residue integer index for this specific atom
            res_idx = res_indices[i]
            
            # 2. Extract the matching residue string name (fallback to 'MOL' if missing)
            res_name = res_names_dict.get(res_idx, "MOL")
            
            mol.add_atom(a,resname=res_name, resnum=res_idx+1)
            
        for x in bonds:
            mol.bonds.append(pmd.Bond(alist[x[0]],alist[x[1]]))
            
        return mol

    
    def GetASEAtoms(self):
        """
        Generate an ASE Atoms instance of the molecule.

        Returns
        -------
        atoms : ase.Atoms
            A newly initialized ASE structural node tracking properties.
        """
        import numpy as np
        import ase
        
        crds = np.array(self.data["positions"])
        eles = self.data["elements"]
        #q    = self.data["charge"]
        q    = self.GetCharge()
        m    = self.data["spin"]
        qs   = self.data["charges"]
        atlist = "".join( ["%s1"%(ele) for ele in eles ] )

        atoms = ase.Atoms(atlist,positions=crds,charges=qs)
        atoms.info["charge"] = q
        atoms.info["spin"] = m
        atoms.calc = None

        return atoms



    def GetRDKitAtoms(self):
        """
        Generate an RDKit Mol instance containing full residue and bonding details
        by locking spatial stereochemistry before the bond-order optimization pass.
        """
        import numpy as np
        from rdkit import Chem
        from rdkit import rdBase
        from rdkit.Chem import rdDetermineBonds
        from rdkit.Chem import AtomPDBResidueInfo

        eles = self.data["elements"]
        bonds = self.data["bonds"]
        crds = np.array(self.data["positions"])
        
        rdBase.DisableLog('rdApp.warning')
        rdBase.DisableLog('rdApp.error')

        # 1. Create a blank template tracking your strict 0-based sequence
        mol = Chem.Mol()
        edit_mol = Chem.EditableMol(mol)
        for symbol in eles:
            r_atom = Chem.Atom(symbol)
            r_atom.SetFormalCharge(0)
            edit_mol.AddAtom(r_atom)
        
        # 2. Inject your bare topology index pairs as explicit single bonds first.
        # This completely guarantees the core ring sizes stay true to your input file definitions.
        for idx1, idx2 in bonds:
            edit_mol.AddBond(idx1, idx2, Chem.BondType.SINGLE)
        mol = edit_mol.GetMol()

        # 3. Inject the true 3D Coordinates using your numpy array directly
        conf = Chem.Conformer(len(eles))
        for i in range(len(eles)):
            conf.SetAtomPosition(i, (float(crds[i, 0]), float(crds[i, 1]), float(crds[i, 2])))
        mol.AddConformer(conf)

        # 4. Perceive connectivity and lock spatial stereochemistry FIRST.
        # By assigning stereochemistry on the locked single-bonded framework, RDKit captures
        # the true coordinate priorities before the bond-order optimizer alters the valences.
        rdDetermineBonds.DetermineConnectivity(mol)
        Chem.AssignStereochemistryFrom3D(mol)

        # 5. Native Bond Type Optimization Engine
        # With the spatial stereocenters securely embedded, the valence engine resolves
        # carbonyls, sulfones, and aromatic paths without mutating the adjacent priorities.
        try:
            rdDetermineBonds.DetermineBondOrders(
                mol, 
                charge=self.GetCharge(), 
                allowChargedFragments=True,
                embedChiral=False
            )
        except ValueError:
            mol.UpdatePropertyCache(strict=False)

        # 6. Standard structural validation and cleanup pass
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)

        # --- Inject Residue Information into the final RDKit Molecule ---
        # FIXED: Cleared the faulty asterisk fallback syntax completely.
        res_indices = self.data.get("residxs", [0] * len(eles))
        res_names_dict = self.data.get("resnames", {0: "MOL"})
        names = self.data.get("names", eles)

        for i, atom in enumerate(mol.GetAtoms()):
            if i < len(eles):
                atom_name = names[i]
                res_idx = int(res_indices[i])
                res_name = res_names_dict.get(res_idx, "MOL")
                res_num = res_idx + 1

                res_info = AtomPDBResidueInfo(
                    atomName=atom_name.ljust(4),
                    residueName=res_name[:3].upper(),
                    residueNumber=res_num,
                    chainId="A"
                )
                atom.SetMonomerInfo(res_info)

        rdBase.EnableLog('rdApp.warning')
        rdBase.EnableLog('rdApp.error')

        return mol


    def Update(self,ene,crds,frcs):
        """
        Update the molecule's potential energy, coordinates, and forces.

        Parameters
        ----------
        ene : float
            The updated molecular potential energy valued in eV.
        crds : list of list of float or numpy.ndarray
            An (N, 3) geometry matrix tracking updated positions in Angstroms.
        frcs : list of list of float or numpy.ndarray or None
            An (N, 3) matrix mapping atomic forces in eV/Angstrom.

        Returns
        -------
        None
        """
        self.data["energy"] = ene
        try:
            self.data["positions"] = crds.tolist()
        except Exception:
            self.data["positions"] = crds
        self.data["forces"] = None
        if frcs is not None:
            try:
                self.data["forces"] = frcs.tolist()
            except Exception:
                self.data["forces"] = frcs


    def clone_geometry(self, coords=None, ene=None, frcs=None):
        """Return a topology-sharing clone with independent ``data`` / positions.

        Shares immutable-ish topology lists (``elements``, ``bonds``, ...) via a
        shallow ``data`` dict copy. Clears the NetworkX ``graph`` cache. Use this
        instead of ``deepcopy`` when only energy / coordinates / forces change
        (wavefront worker results).
        """
        import copy

        out = copy.copy(self)
        out.data = dict(self.data)
        if hasattr(out, "graph"):
            out.graph = None
        use_ene = out.data.get("energy") if ene is None else ene
        use_crds = out.data.get("positions") if coords is None else coords
        use_frc = out.data.get("forces") if frcs is None else frcs
        if coords is not None or ene is not None or frcs is not None:
            out.Update(use_ene, use_crds, use_frc)
        return out


    def SaveCrds(self,fname_or_fh,fmt=None):
        """
        Export the structural coordinates and topology details to disk.

        Parameters
        ----------
        fname_or_fh : str or file-like object
            The disk filepath string destination or an opened file descriptor handle.
        fmt : str, default None
            Explicit file format designation ('XYZ', 'MOL2', 'RST7', 'PDB'). If 
            None, the target format family is determined from file suffixes.

        Returns
        -------
        None

        Raises
        ------
        Exception
            If an unexpected handle context is passed without an explicit format flag,
            or if an unknown file suffix is encountered.
        """
        
        from pathlib import Path
        import ase.io
        import numpy as np

        import parmed as pmd
        from parmed.structure import Structure
        from parmed.topologyobjects import Atom, Bond, Residue
        from . constants import GetAtomicNumber
        
        eles  = self.data["elements"]
        types = self.data["types"]
        qs    = self.data["charges"]
        #q     = self.data["charge"]
        q     = self.GetCharge()
        crds  = self.data["positions"]
        bonds = self.data["bonds"]

        if hasattr(fname_or_fh, 'read') or hasattr(fname_or_fh, 'write'):
            if fmt is None:
                raise Exception("Expected fmt because fname_or_fh is a file descriptor")
            fh = fname_or_fh
        else:
            suffix = Path(fname_or_fh).suffix
            if suffix == ".xyz":
                fmt="XYZ"
                fh = fname_or_fh
            elif suffix == ".mol2":
                fmt="MOL2"
                fh = fname_or_fh
            elif suffix == ".rst7":
                fmt="RST7"
                fh = open(fname_or_fh,"w")
            elif suffix == ".pdb":
                fmt="PDB"
                fh = open(fname_or_fh,"w")
            elif suffix == ".json":
                fmt="JSON"
                fh = open(fname_or_fh,"w")
            else:
                raise Exception(f"Unknown suffix {suffix}")
        fmt = fmt.upper()
        if fmt == "XYZ":
            atoms = self.GetASEAtoms()
            ase.io.write(fh,atoms,format="extxyz")
        elif fmt == "MOL2":
            mol = self.GetParmedAtoms()
            mol.save(fh,overwrite=True,format="MOL2")
        elif fmt == "RST7":
            crds = np.array(crds)
            #fh = open(fname,"w")
            write_amber_restart(crds, fh)
            #fh.close()
        elif fmt == "PDB":
            mol = self.GetParmedAtoms()
            mol.save(fh,overwrite=True,format="PDB")
        elif fmt == "JSON":
            import json
            #los = ListOfStruct( [self] )
            if self.data["constraints"] is None:
                self.data["constraints"] = []
            json.dump([ self.data ],fh,indent=4)
        else:
            raise Exception(f"Unknown format {fmt}")

        
    def GetInchi(self,tautdep=True):
        return RDKitMol2Inchi(self.GetRDKitAtoms(),tautdep=tautdep)

    def GetInchiToken(self,tautdep=True):
        return RDKitMol2InchiToken(self.GetRDKitAtoms(),tautdep=tautdep)

    def GetInchiKey(self,tautdep=True):
        return RDKitMol2InchiKey(self.GetRDKitAtoms(),tautdep=tautdep)

    def GetSmiles(self,tautdep=True):
        return RDKitMol2Smiles(self.GetRDKitAtoms(),tautdep=tautdep)

    def GetIUPAC(self):
        from nispo import smiles_to_iupac
        mol = self.GetSmiles(tautdep=True)
        with silence_stderr():
            s = smiles_to_iupac( mol )
        return s
    
    def GetGraph(self):
        """
        Generate a search graph mapping covalent bonding connectivity networks.

        Returns
        -------
        graph_copy : ffpopt.GraphSearch
            A deep copy of the molecular bond connectivity search graph.
        """
        
        from . AmberParm import bonds2graph
        from copy import deepcopy
        if self.graph is None:
            self.graph = bonds2graph(self.data["bonds"])
        return deepcopy(self.graph)

    
    def GetFunctionalGroups(self):
        """
        Detect and catalog localized chemical functional groups.

        Returns
        -------
        fg_dict : dict
            A dictionary mapping functional group classification strings (keys)
            to arrays of atomic indices tracking their boundaries (values).
        """
        from . FindFuncGrps import FindFuncGrps_from_rdkit
        return FindFuncGrps_from_rdkit( self.GetRDKitAtoms() )

    
    def GetFunctionalGroupCharges(self):
        """
        Compute the net partial charge sum of detected functional groups.

        Returns
        -------
        fgq : dict
            A dictionary mapping functional group names (keys) to sequences
            of floating-point values tracking their net sums (values).
        """
        from collections import defaultdict as ddict
        fg = self.GetFunctionalGroups()
        fgq = ddict(list)
        for key in fg:
            fgq[key] = [ self.GetChargeSum(idxs)
                         for idxs in fg[key] ]
        return fgq

    
    def GetChargeSum(self,idxs):
        """
        Compute the partial charge sum across a specific set of atom indices.

        Parameters
        ----------
        idxs : list of int
            An array of 0-based atomic row indices to subset.

        Returns
        -------
        sumqs : float
            The net partial sum of the selected atomic charges.
        """
        qs = [ self.data["charges"][idx] for idx in idxs ]
        return sum(qs)


    def GetHilfikerCharges(self):
        """
        Compute atomic charges using the Hilfiker fitting pipeline scheme.

        Returns
        -------
        charges : list of float
            An ordered array of calculated atomic partial charges.
        """
        from . RespFit import hilfiker_charges
        atoms = self.GetASEAtoms()
        q = int(round(sum(self.data["charges"])))
        return hilfiker_charges(atoms,total_charge=q)

    
    def GetEspalomaCharges(self):
        """
        Compute atomic partial charges using the Espaloma deep learning engine.

        Returns
        -------
        charges : list of float
            An ordered sequence of machine-learning predicted partial charges.
        """
        import sys
        import warnings
        from types import ModuleType
        from functools import partial
        import torch
        import rdkit.rdBase as rdb
        from . RespFit import espaloma_charge

        # 1. Suppress Python-level warnings
        warnings.filterwarnings("ignore", category=UserWarning, module="torchdata")
        warnings.filterwarnings("ignore", message=".*Recommend creating graphs.*")

        # 2. Inject dummy graphbolt structure into cache if not present
        if 'dgl.graphbolt' not in sys.modules:
            fake_graphbolt = ModuleType('dgl.graphbolt')
            fake_graphbolt.__path__ = []
            sys.modules['dgl.graphbolt'] = fake_graphbolt

        # 3. Override PyTorch pickler settings for legacy model unpickling
        torch.load = partial(torch.load, weights_only=False)

        # 4. Safely import your frameworks inside the function context
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        from espaloma_charge import charge
        
        molecule = self.GetRDKitAtoms()
        #print("Total charge:", Chem.GetFormalCharge(molecule))
        #print("Expected charge:",int(round(self.data["charge"])))
        
        # 5. Temporarily silence RDKit C++ streams only during the charge call
        # This protects the rest of your ffpopt pipeline's error visibility
        rdb.DisableLog('rdApp.warning')
        rdb.DisableLog('rdApp.error')
        rdb.DisableLog('rdApp.info')

        try:
            charges = espaloma_charge(molecule,total_charge=self.GetCharge())
        finally:
            # 6. Re-enable the C++ log streams before control exits the method
            rdb.EnableLog('rdApp.warning')
            rdb.EnableLog('rdApp.error')
            rdb.EnableLog('rdApp.info')

        #print("sum=",sum(charges))
            
        return charges

    
    def set_dihedral(self,a,b,c,d,value):
        """
        Modify a specific dihedral angle and update atomic positions.

        Parameters
        ----------
        a : int
            The 0-based atomic row index tracking the first node of the dihedral sequence.
        b : int
            The 0-based atomic row index tracking the second node of the dihedral sequence.
        c : int
            The 0-based atomic row index tracking the third node of the dihedral sequence.
        d : int
            The 0-based atomic row index tracking the fourth node of the dihedral sequence.
        value : float
            The target geometric value assigned to the angle in degrees.

        Returns
        -------
        None
        """
        from . Constraints import Constraint
        atoms = self.GetASEAtoms()
        con = Constraint("dihed",[a,b,c,d],
                         value=value,
                         graph=self.GetGraph())
        atoms = con.modify(atoms)
        self.Update(0,atoms.get_positions(),None)
 
        
    
    #
    # These mimic methods of an ASE Atoms object
    #
    
    def get_dihedral(self,a,b,c,d):
        """
        Calculate the dihedral angle defined by four specific atom indices.

        Parameters
        ----------
        a : int
            The first 0-based atom index.
        b : int
            The second 0-based atom index.
        c : int
            The third 0-based atom index.
        d : int
            The fourth 0-based atom index.

        Returns
        -------
        angle : float
            The calculated torsion value expressed in degrees.
        """
        return self.GetASEAtoms().get_dihedral(a,b,c,d)

    def get_angle(self,a,b,c):
        """
        Calculate the valence angle defined by three specific atom indices.

        Parameters
        ----------
        a : int
            The first 0-based atom index.
        b : int
            The vertex 0-based atom index.
        c : int
            The third 0-based atom index.

        Returns
        -------
        angle : float
            The calculated angle value expressed in degrees.
        """
        return self.GetASEAtoms().get_angle(a,b,c)

    def get_distance(self,a,b):
        """
        Calculate the Euclidean bond length between two specific atom indices.

        Parameters
        ----------
        a : int
            The first 0-based atom index.
        b : int
            The second 0-based atom index.

        Returns
        -------
        bond_length : float
            The physical distance value expressed in Angstroms.
        """
        return self.GetASEAtoms().get_distance(a,b)

    def get_potential_energy(self):
        """
        Retrieve the stored potential energy of the molecular layout.

        Returns
        -------
        potential_energy : float or None
            The recorded structural layout energy valued in eV.
        """
        return self.data["energy"]

    def get_forces(self):
        """
        Retrieve the stored atomic forces matrix.

        Returns
        -------
        forces : numpy.ndarray
            An (N, 3) matrix mapping precise atomic force arrays in eV/Angstrom.
        """
        import numpy as np
        return np.array(self.data["forces"])

    def get_positions(self):
        """
        Retrieve a copy of the coordinate matrix.

        Returns
        -------
        coords : numpy.ndarray
            An (N, 3) matrix mapping precise Cartesian coordinates in Angstroms.
        """
        import numpy as np
        return np.array(self.data["positions"])

    def copy(self):
        """
        Generate a deep copy of this Struct instance.

        Returns
        -------
        struct_copy : Struct
            An independent, completely decoupled duplicate of the Struct object.
        """
        import copy
        return copy.deepcopy(self)

    
class ListOfStruct(object):
    """
    Collection manager and orchestrator for multiple Struct instances.

    This class serves as a high-level container for a list of `Struct` objects,
    enabling bulk operations such as serializing/deserializing batches of molecular
    conformations, searching structures by identification tags, and orchestrating
    ASE calculator assignments across structural frames.

    Attributes
    ----------
    SANDERMODES : list of str
        Supported execution string modes identifying valid AMBER sander force field or
        semi-empirical backends ('SANDER', 'DFTB3', 'DFTB2', 'AM1D').
    structs : list of Struct
        The internal list containing the wrapped Struct instances.
    args : dict or object or None
        Active command line arguments or dictionary configuration blocks used to determine
        calculator variables and environmental settings.
    calc : ase.calculators.calculator.Calculator or None
        The currently cached or assigned backend evaluation calculator reference.
    """
    
    
    SANDERMODES=["SANDER","DFTB3","DFTB2","AM1D"]
    
    def __init__(self,structs):
        """
        Initialize a ListOfStruct collection wrapper instance.

        Parameters
        ----------
        structs : list of Struct
            A list containing the Struct objects to pack inside this collection.

        Raises
        ------
        TypeError
            If the supplied structs parameter is a string instead of an iterable collection.
        """
        from copy import deepcopy

        if isinstance(structs, str):
            raise TypeError(f"Expected a list or dict, but received a string: {structs!r}")
        
        self.structs = deepcopy(structs)
        self.args  = None
        self.calc  = None

        
    def __getstate__(self):
        """
        Extract the serialization state dictionary while omitting active calculator handles.

        Returns
        -------
        state : dict
            The object's instance state dictionary with the 'calc' reference set to None.
        """
        import copy
        # Copy the object's state dictionary
        state = self.__dict__.copy()
        # Remove the specific attribute
        if 'calc' in state:
            #del state['calc']
            state['calc'] = None
        return state
        
    @classmethod
    def from_file(cls,filename):
        """
        Instantiate a ListOfStruct instance from a JSON dataset or a Tripos MOL2 structural file.

        Parameters
        ----------
        filename : str
            The absolute or relative file path targeting the input JSON array or MOL2 file.

        Returns
        -------
        ListOfStruct
            A newly initialized ListOfStruct collection containing the parsed structures.

        Raises
        ------
        Exception
            If both JSON decoding and MOL2 file parsing fail to read the target file.
        """

        ss = None
        try:
            import json
            fh = open(filename,'r')
            data = json.load(fh)
            ss = [Struct(x) for x in data]
        except Exception as JsonError:
            try:
                ss = [ Struct.from_mol2(filename) ]
            except Exception as Mol2Error:
                try:
                    prefix,charge = parse_filename_charge(filename)
                    ss = [ Struct.from_xyz(prefix,charge) ]
                except Exception as XYZError:
                    try:
                        mol = GetRDKitMolFromString(filename,quiet=True)
                        ss = [ Struct.from_rdkit(mol) ]
                    except Exception as SmilesError:
                        raise Exception(JsonError)
        return cls( ss )

    

    def __iter__(self):
        """
        Return an iterator over the packed Struct collection.

        Returns
        -------
        iterator
            An iterator traversing the internal list of packed molecular frames.
        """
        return iter(self.structs)

    
    def __len__(self):
        """
        Retrieve the total number of Struct instances packed inside this collection.

        Returns
        -------
        length : int
            The total index count of packed structural configuration frames.
        """
        return len(self.structs)

    
    def __setitem__(self, index, value):
        """
        Assign or overwrite a Struct instance at a specific index boundary location.

        Parameters
        ----------
        index : int or slice
            The target placement index marker position to update.
        value : Struct
            The replacement Struct object template to bind onto the index position.
        """
        self.structs[index] = value

        
    def __getitem__(self, index):
        """
        Extract a specific Struct conformation frame by index.

        Parameters
        ----------
        index : int or slice
            The index marker tracking the targeted structural record block.

        Returns
        -------
        struct : Struct
            The specific molecular conformation instance localized at the index path.
        """
        # This maps indexing directly to the list member
        return self.structs[index]

    
    def save(self,filename):
        """
        Serialize and dump the entire list of packed structural records to disk as a JSON matrix.

        Parameters
        ----------
        filename : str
            The destination target path string where the packed JSON array is written.
        """
        import json


        if filename.endswith(".json"):
            fh = open(filename,"w")
            for x in self.structs:
                if x.data["constraints"] is None:
                    x.data["constraints"] = []
            json.dump([ x.data for x in self.structs ],fh,indent=4)
        else:
            from pathlib import Path
            names = []
            opath = Path(filename)
            for inp in self.structs:
                names.append(inp.data["name"])
                if len(self.structs) > 1:
                    base = opath.with_suffix('')
                    s = opath.suffix
                    s = "_" + inp.data["name"] + str(s)
                    o = str(base) + s
                else:
                    o = str(opath)
                inp.SaveCrds(o)



    def GetByName(self,name):
        """
        Search the collection and find a Struct instance matching a unique text label name.

        Parameters
        ----------
        name : str
            The alphanumeric identity name string tracking the requested molecule.

        Returns
        -------
        matched_struct : Struct or None
            The targeted Struct instance object matching the lookup criteria, or None if missing.
        """
        keep = None
        for s in self.structs:
            if s.data["name"] == name:
                keep=s
                break
        return keep
        

        
    def SetArgs(self,args):
        """
        Configure system properties, device execution targets, and external tool keyword parameters.

        Parameters
        ----------
        args : dict or object
            An argument tracking context or configuration dictionary containing setup metrics. 
            Supported fields include 'cpu', 'geometric_coordsys', and 'psi4_memory'.
        """
        from copy import deepcopy
        self.args = deepcopy(args)

        if hasattr(self.args,"cpu"):
            if self.args.cpu:
                import os
                os.environ['JAX_PLATFORMS'] = 'cpu'
                os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

        geometric_kwargs = {}
        psi4_kwargs = {}

        if isinstance(args,dict):
            for key in ["coordsys","maxiter","converge","enforce"]:
                if key in args:
                    geometric_kwargs[key] = args[key]
            for key in ["memory","num_threads"]:
                if key in args:
                    geometric_kwargs[key] = args[key]
    
        else:
            try:
                geometric_kwargs = { "coordsys": str(args.geometric_coordsys),
                                     "maxiter": str(args.geometric_maxiter),
                                     "converge": str(args.geometric_converge),
                                     "enforce": str(args.geometric_enforce) }
            except:
                pass

            try:
                psi4_kwargs = { "memory": str(args.psi4_memory),
                                "num_threads": str(args.psi4_num_threads) }
            except:
                pass

        
        self.extra_args = { "geometric": geometric_kwargs,
                            "psi4": psi4_kwargs }

        #self.calc = self.BuildCalc(self.structs[0])
        self.calc = None

        
            
    def BuildCalc(self,struct,calc=None):
        """
        Build or recycle an active potential evaluator matching the target frame's settings.

        Parameters
        ----------
        struct : Struct
            The molecular conformation layout whose parameters dictate the calculator's setup.
        calc : ase.calculators.calculator.Calculator, default None
            An optional preexisting calculator context to reuse if topological states are matching.

        Returns
        -------
        assigned_calculator : ase.calculators.calculator.Calculator
            The newly compiled or recycled evaluation calculator assigned to the structure.
        """
        from . ase.calculator import GenCalculator

        if calc is None:
            if self.calc is not None:
                if isinstance(self.args, dict):
                    want_mode = str(self.args.get("model", "SANDER")).upper()
                else:
                    want_mode = str(getattr(self.args, "model", "SANDER")).upper()
                if self.calc.parm == struct.data["parm"] and \
                   self.calc.spin == struct.data["spin"] and \
                   self.calc.charge == struct.GetCharge() and \
                   str(getattr(self.calc, "mode", "")).upper() == want_mode:
                    calc = self.calc
                    #print("CALC: reuse previous")

        model = None
        if calc is None:
            if isinstance(self.args,dict):
                if "psi4_memory" in self.args:
                    memory = self.args["psi4_memory"]
                elif "memory" in self.args:
                    memory = self.args["memory"]
                else:
                    memory = "1gb"
                if "psi4_num_threads" in self.args:
                    num_threads = int(self.args["psi4_num_threads"])
                elif "num_threads" in self.args:
                    num_threads = int(self.args["num_threads"])
                else:
                    num_threads = 1
                    
                if "mfile" in self.args:
                    mfile = self.args["mfile"]
                else:
                    mfile = None
                if "model" in self.args:
                    model = self.args["model"]
            else:
                memory = getattr(self.args,"psi4_memory","1gb")
                num_threads = getattr(self.args,"psi4_num_threads",1)
                mfile = None
                if hasattr(self.args,"mfile"):
                    mfile = self.args.mfile
                model = "SANDER"
                if hasattr(self.args,"model"):
                    model = self.args.model
                
            kwargs = { "memory": memory,
                       "num_threads": num_threads,
                       "mfile": mfile }
            
            charge = struct.GetCharge()
            spin = struct.data["spin"]
            parm = struct.data["parm"]
            mode = model.upper()
            
            if mode in ListOfStruct.SANDERMODES:
                import tempfile
                import os
                import numpy as np
                from pathlib import Path
                
                tmpfile_loc = "./tmpfiles"
                if not Path(tmpfile_loc).is_dir():
                    os.makedirs(tmpfile_loc, exist_ok=True)

                fd, path = tempfile.mkstemp(dir=tmpfile_loc,suffix=".rst7")
                fh = os.fdopen(fd,'w')
                write_amber_restart( np.array(struct.data["positions"]), fh )
                fh.close()

                try:
                    import sander
                    sander.cleanup()
                except:
                    pass

                calc = GenCalculator(mode,charge,spin,parm,path,**kwargs)
                os.remove(path)
            else:
                calc = GenCalculator(mode,charge,spin,parm,None,**kwargs)
        self.calc = calc
        return calc

    
    def BuildRestrainedCalc(self,s,reslist=None):
        """
        Build an ASE calculator coupled to geometric restraints.

        Parameters
        ----------
        s : Struct
            The target conformation instance whose variables map to the base parameters.
        reslist : RestraintList, default None
            Optional explicit sequence of restraints. If None, it automatically falls back
            to pulling structural records stored directly inside `s.restraints`.

        Returns
        -------
        restrained_calculator : ase.calculators.calculator.Calculator
            The active RestrainedCalculator instance object wrapping the target potential functions.
        """
        calc = self.BuildCalc(s)
        rcalc = calc

        res = None
        if reslist is not None:
            res = reslist
        elif s.restraints is not None:
            res = s.restraints
            
        if res is not None:
            from . ase.calculator import RestrainedCalculator
            rcalc = RestrainedCalculator(calc,res)
            
        return rcalc



from . ase.calculator import GenCalculator
from ase.calculators.calculator import Calculator,all_changes

class RestCalculator(GenCalculator):
    
    implemented_properties = ['energy','forces','free_energy']
    nolabel=True
    def __init__(self,mode,inp,**kwargs):
        los = ListOfStruct.from_file(inp)
        kwargs["model"] = mode.upper()
        #print("kwargs=",kwargs)
        los.SetArgs(kwargs)
        self.charge = los.structs[0].GetCharge()
        #print("kwargs=",kwargs)
        self.calc = los.BuildRestrainedCalc(los.structs[0])
        #print("kwargs=",kwargs)
        Calculator.__init__(self,**kwargs)

    def calculate(self,
                  atoms=None,
                  properties=None,
                  system_changes=all_changes):
        import numpy as np
        import ase
        if self.charge is not None:
            atoms.info["charge"] = self.charge
        if properties is None:
            properties = self.implemented_properties
        Calculator.calculate(self, atoms, properties, system_changes)
        atoms.calc =  self.calc
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces




