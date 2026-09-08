from ase.calculators.calculator import Calculator, all_changes

class RDKitUFFCalculator(Calculator):
    """
    An ASE Calculator that dynamically manages an RDKit UFF engine.

    This calculator builds an RDKit RWMol representation and a Universal
    Force Field (UFF) engine on the fly from the incoming ASE Atoms object.
    It preserves the underlying chemical topology to minimize overhead, while
    safely rebuilding the force field evaluator at each step to prevent coordinate
    disconnection bugs during geometry optimizations.

    Parameters
    ----------
    charge : int, default 0
        The total net electronic charge of the molecular system. This value is
        passed directly to the RDKit bond-order perception engine to resolve
        correct valences and formal charges.
    **kwargs : dict
        Additional keyword arguments passed straight to the base ASE
        `Calculator` constructor.

    Attributes
    ----------
    charge : int
        The total net charge of the molecular system.
    uff : rdkit.ForceField.rdForceField.ForceField or None
        The active RDKit UFF force field instance, or None if not initialized.
    mol : rdkit.Chem.rdchem.Mol or None
        The persistent reference to the active RDKit molecule structure.
    implemented_properties : list of str
        Properties that this calculator is capable of computing:
        `['energy', 'forces']`.

    See Also
    --------
    ase.calculators.calculator.Calculator : The base ASE calculator class.
    """
    implemented_properties = ['energy', 'forces']

    def __init__(self, charge=0, **kwargs):
        super().__init__(**kwargs)
        self.charge = int(charge)
        self.uff = None
        self.mol = None  # Persistent storage for the active RDKit molecular graph
        self._chemical_symbols = []  # Internal tracking for structural changes
        self._previous_forces = None  # Diagnostic memory cache to test force variation

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """
        Compute the UFF potential energy and atomic forces.

        Synchronizes the atomic positions from the ASE Atoms object to the
        underlying RDKit conformation, rebuilds the force field evaluator to capture
        the updated geometry, and runs the calculation. If the number of atoms or
        sequence of elements changes relative to a previous call, the molecular graph
        is cleanly re-perceived on the fly.

        Parameters
        ----------
        atoms : ase.Atoms or None, default None
            The ASE Atoms object to evaluate. If None, the previously stored
            Atoms object is evaluated.
        properties : list of str or None, default None
            The specific properties requested for calculation. Defaults to
            `['energy', 'forces']`.
        system_changes : list of str, default all_changes
            List of structural properties that have changed since the last
            evaluation, used by the base class.
        """
        if properties is None: 
            properties = ['energy', 'forces']
        super().calculate(atoms, properties, system_changes)

        import numpy as np
        from rdkit import Chem
        from rdkit.Geometry import Point3D

        current_symbols = self.atoms.get_chemical_symbols()
        num_atoms = len(self.atoms)
        positions = self.atoms.get_positions()

        # Rebuild the core chemical topology ONLY if elements or size change
        if self.mol is None or self._chemical_symbols != current_symbols:
            self.uff = None 
            self.mol = None
            self._previous_forces = None
            self._chemical_symbols = list(current_symbols)

            # Create a mutable RWMol object as expected by RDKit's C++ signatures
            rw_mol = Chem.RWMol()
            for symbol in current_symbols:
                atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(symbol)
                atom = Chem.Atom(atomic_num)
                atom.SetNoImplicit(True)
                rw_mol.AddAtom(atom)
            
            # Setup a basic spatial conformer template matching the RWMol size
            conformer = Chem.Conformer(num_atoms)
            conformer.Set3D(True)
            
            # Populate with the initial positions to accurately perceive bonds
            for i, (x, y, z) in enumerate(positions):
                conformer.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
            rw_mol.AddConformer(conformer, assignId=True)

            # Perceive connectivity and bond orders dynamically on the RWMol
            from rdkit.Chem import rdDetermineBonds
            try:
                rdDetermineBonds.DetermineBonds(rw_mol, charge=int(self.charge))
            except Exception as e:
                raise RuntimeError(f"RDKit DetermineBonds execution failed: {str(e)}")

            # Cache the generated molecular graph permanently
            self.mol = rw_mol
        # 2. Push updated dynamic positions directly into our persistent molecule object
        from rdkit.Geometry import Point3D
        import sys
        conf = self.mol.GetConformer()
        for i, (x, y, z) in enumerate(positions):
            conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

        # 3. Force field lookup instantiation block matching active pointer address layouts
        from rdkit.Chem import AllChem
        self.uff = AllChem.UFFGetMoleculeForceField(self.mol)
        if self.uff is None:
            raise ValueError(
                f"Failed to initialize RDKit UFF for the provided system "
                f"with net charge {self.charge}."
            )
        self.uff.Initialize()

        # 4. Handle conversion factors cleanly using native ASE utilities
        from ase.units import kcal, mol as ase_mol
        kcal_mol_to_ev = kcal / ase_mol

        # 5. Compute potential energy and gradients
        import numpy as np
        self.results['energy'] = self.uff.CalcEnergy() * kcal_mol_to_ev
        
        grad_kcal = np.array(self.uff.CalcGrad())
        forces_ev = -grad_kcal.reshape((-1, 3)) * kcal_mol_to_ev

        # 6. SILENT RUNTIME HYPOTHESIS INTEGRITY CHECK
        # Keeps calculation safe without cluttering your optimization stdout streams
        if self._previous_forces is not None:
            force_diff = np.max(np.abs(forces_ev - self._previous_forces))
            
            if force_diff < 1e-6:
                tracker_url = ".".join(["github", "com"]) + "/".join(["", "rdkit", "rdkit", "issues"])
                report = [
                    "\n" + "="*80,
                    "CRITICAL ERROR: RUNTIME STATIC FORCE RUNAWAY DETECTED",
                    "="*80,
                    f"Max Force Delta : {force_diff:.8f} eV/Å (Less than 1e-6 tolerance bounds)",
                    f"Action Required : Verify model parameters or coordinate trackers via {tracker_url}",
                    "="*80 + "\n"
                ]
                print("\n".join(report))
                sys.stdout.flush()
                raise RuntimeError("Static force loop anomaly caught during optimization runtime execution.")

        # Cache the current step forces into memory for the next structural assertion test
        self._previous_forces = np.copy(forces_ev)
        self.results['forces'] = forces_ev
