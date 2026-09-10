#!/usr/bin/env python3



def DeleteDihedrals(p,list_of_idxs):
    """
    Delete dihedrals from the Parm object p.
    
    Parameters
    ----------
    p : parmed.AmberParm
        The Parm object from which dihedrals will be deleted.
    list_of_idxs : list of tuples
        A list of tuples, where each tuple contains 4 integers representing the indices of the atoms forming the dihedral to be deleted.
    
    Returns
    -------
    None
    
    """
    from parmed.tools.actions import deleteDihedral, addDihedral
    for idxs in list_of_idxs:
        cmd = deleteDihedral(p,
                             f"@{idxs[0]+1}",
                             f"@{idxs[1]+1}",
                             f"@{idxs[2]+1}",
                             f"@{idxs[3]+1}")
    
        cmd.execute()
        
    
    

def ChangeDihedrals(p,idxs,xs,fc=None,bytype=False):
    """ Change dihedrals in the Parm object p.

    Parameters
    ----------
    p : parmed.AmberParm
        The Parm object in which dihedrals will be changed.
    idxs : list of int
        A list of 4 integers representing the indices of the atoms forming the dihedral to be changed.
    xs : list of DihedralType
        A list of DihedralType objects representing the new dihedral parameters.
    fc : float, optional
        A float value to set the force constant for all dihedrals. If None, the force constants from xs will be used.
    bytype : bool, optional
        If True, the function will change dihedrals by type instead of by indices. Default is False.
    
    Returns
    -------
    None

    """

    from collections import defaultdict as ddict
    from parmed.tools.actions import deleteDihedral, addDihedral

    if bytype:

        ftypes = tuple([p.atoms[idx].type for idx in idxs])
        rtypes = tuple(list(ftypes)[::-1])
        allidxs = []
        for x in p.dihedrals:
            kidxs = (x.atom1.idx,x.atom2.idx,x.atom3.idx,x.atom4.idx)
            tidxs = (x.atom1.type,x.atom2.type,x.atom3.type,x.atom4.type)
            if tidxs == ftypes:
                allidxs.append(kidxs)
            elif tidxs == rtypes:
                allidxs.append( tuple(list(kidxs)[::-1]) )
        allseles = [ [ f"@{idx[i]+1}" for i in range(4) ]
                     for idx in allidxs ]

    else:
        allseles = [ [ f"@{i+1}" for i in idxs ] ]

        
    for seles in allseles:
        
        cmd = deleteDihedral(p,seles[0],seles[1],seles[2],seles[3])    
        cmd.execute()

        for x in xs:
            per = x.per
            phase = x.phase
            scee = x.scee
            scnb = x.scnb
            
            k = x.phi_k
            if fc is not None:
                k = fc
    
            cmd = addDihedral\
                (p,
                 seles[0],seles[1],seles[2],seles[3],
                 k,per,phase,scee,scnb)
            cmd.execute()

    

    

class PrimDihedFcn(object):
    """ A class representing a single dihedral function with a force constant, phase, and period.
    
    Parameters
    ----------
    fc : float
        The force constant for the dihedral function.
    phase : float   
        The phase of the dihedral function in degrees.
    per : int
        The period of the dihedral function.
    
    Attributes
    ----------
    fc : float
        The force constant for the dihedral function.
    phase : float
        The phase of the dihedral function in degrees.
    per : int
        The period of the dihedral function.
    
    """
    def __init__(self,fc,phase,per):
        self.fc = fc
        self.phase = phase
        self.per = per

    def __str__(self):
        """ Return a string representation of the PrimDihedFcn object. """
        return f"({self.fc},{self.phase},{self.per})"
        
    def CptEne(self,ang):
        """ Calculate the energy contribution of the dihedral function for a given angle.
        
        Parameters
        ----------
        ang : float
            The angle in degrees for which the energy contribution is calculated.

        Returns
        -------
        float
            The energy contribution of the dihedral function for the given angle.
        
        """
        return self.fc * self.CptEterm(ang)

    def CptEterm(self,ang):
        """ Calculate the energy term of the dihedral function for a given angle.

        Parameters
        ----------
        ang : float
            The angle in degrees for which the energy term is calculated.
        
        Returns
        -------
        float
            The energy term of the dihedral function for the given angle.
        
        """
        import numpy as np
        a = (self.per * ang + self.phase)*(np.pi/180)
        return 1+np.cos(a)

    
class MultiDihedFcn(object):
    """ A class representing a multi-term dihedral function.
    
    Parameters
    ----------
    idxs : list of int
        A list of 4 integers representing the indices of the atoms forming the dihedral.
    prims : list of PrimDihedFcn
        A list of PrimDihedFcn objects representing the individual terms of the dihedral function.
    
    Attributes
    ----------
    idxs : list of int
        A list of 4 integers representing the indices of the atoms forming the dihedral.
    prims : list of PrimDihedFcn
        A list of PrimDihedFcn objects representing the individual terms of the dihedral function.


    """
    def __init__(self,idxs,prims):
        import copy
        self.idxs = copy.deepcopy(idxs)
        self.prims = copy.deepcopy(prims)

    def CptEne(self,ang):
        """ Calculate the total energy contribution of the multi-term dihedral function for a given angle.

        Parameters
        ----------
        ang : float
            The angle in degrees for which the total energy contribution is calculated.
        
        Returns
        -------
        float
            The total energy contribution of the multi-term dihedral function for the given angle.
        
        """

        import numpy as np
        outs = []
        for p in self.prims:
            outs.append( p.CptEne(ang) )
        return np.sum(outs,axis=0)

    def SetFCs(self,fcs):
        """ Set the force constants for each term in the multi-term dihedral function.
        
        Parameters
        ----------
        fcs : list of float
            A list of force constants to set for each term in the multi-term dihedral function.
        
        Returns
        -------
        None
        
        Raises
        ------
        Exception
            If the length of fcs does not match the number of terms in the multi-term dihedral function.
            
        """
        if len(fcs) != len(self.prims):
            raise Exception(f"Size mismatch {len(fcs)} {len(self.prims)}")
        for i in range(len(fcs)):
            self.prims[i].fc = fcs[i]


    def SetPhases(self,ps):
        """ Set the force constants for each term in the multi-term dihedral function.
        
        Parameters
        ----------
        ps : list of float
            A list of phases to set for each term in the multi-term dihedral function.
        
        Returns
        -------
        None
        
        Raises
        ------
        Exception
            If the length of ps does not match the number of terms in the multi-term dihedral function.
            
        """
        if len(ps) != len(self.prims):
            raise Exception(f"Size mismatch {len(ps)} {len(self.prims)}")
        for i in range(len(ps)):
            self.prims[i].phase = ps[i]

            
    def ResetFCs(self):
        """ Reset the force constants to default values. """
        self.SetFCs( [1]*len(self.prims) )

    def ResetPhases(self):
        """ Reset the force constants to default values. """
        self.SetPhases( [0]*len(self.prims) )

    def __str__(self):
        """ Return a string representation of the MultiDihedFcn object. """
        return "[%s]"%(",".join([str(x) for x in self.prims]))


    
def CptDihedralEne(atoms,mdfs):
    """ Calculate the dihedral energy for a list of MultiDihedFcn objects.
    
    Parameters
    ----------
    atoms : parmed.amber.AmberParm
        The AmberParm object containing the molecular structure.
    mdfs : list of MultiDihedFcn
        A list of MultiDihedFcn objects representing the dihedral functions.
    
    Returns
    -------
    numpy.ndarray
        A numpy array containing the energy contributions for each MultiDihedFcn object.
    
    """
    import numpy as np
    enes = []
    for mdf in mdfs:
        idxs = mdf.idxs
        ang = atoms.get_dihedral(idxs[0],idxs[1],idxs[2],idxs[3])
        e = mdf.CptEne(ang)
        enes.append(e)
    return np.array(enes)
        



def FindDihedrals(p,idxs,impropers=False):
    """ Find dihedrals in the Parm object p that match the given indices.
    
    Parameters
    ----------
    p : parmed.AmberParm
        The Parm object containing the molecular structure.
    idxs : list of int
        A list of 4 integers representing the indices of the atoms forming the dihedral to be found.
    impropers : bool, optional
        If True, only find improper dihedrals. Default is False.
    
    Returns
    -------
    list of Dihedral
        A list of Dihedral objects that match the given indices.
    
    """
    tkey = tuple(idxs)
    xs=[]
    for x in p.dihedrals:
        if x.improper and not impropers:
            continue
        fkey = (x.atom1.idx,x.atom2.idx,x.atom3.idx,x.atom4.idx)
        rkey = tuple(list(fkey)[::-1])
        #print(tkey,fkey,rkey,fkey == tkey, rkey == tkey, x.atom1.type, x.atom2.type,x.atom3.type,x.atom4.type)
        if fkey == tkey or rkey == tkey:
            xs.append(x)
    return xs




def format_parm_quartet(p, idxs) -> str:
    """Log line: Fourier terms on one 4-atom quartet in ``p``."""

    xs = FindDihedrals(p, idxs)
    if not xs:
        return f"{list(idxs)}: (absent from parm)"
    bits = []
    for x in xs:
        t = x.type
        bits.append(
            f"PK={float(t.phi_k):.4f} n={int(t.per)} γ={float(t.phase):.0f}"
        )
    return f"{list(idxs)}: " + "; ".join(bits)


def summarize_rotors_on_bond(p, idxs):
    """Proper-dihedral terms that share the scanned central bond."""

    a, b = int(idxs[1]), int(idxs[2])
    groups = {}
    for x in p.dihedrals:
        if x.improper:
            continue
        i2, i3 = x.atom2.idx, x.atom3.idx
        if not ((i2 == a and i3 == b) or (i2 == b and i3 == a)):
            continue
        q = (x.atom1.idx, x.atom2.idx, x.atom3.idx, x.atom4.idx)
        key = q if q <= q[::-1] else q[::-1]
        groups.setdefault(key, []).append(x)
    lines = []
    for q, xs in sorted(groups.items()):
        pks = [
            f"n={int(x.type.per)} PK={float(x.type.phi_k):.4f} γ={float(x.type.phase):.0f}"
            for x in xs
        ]
        types = "-".join(p.atoms[i].type for i in q)
        scanned = "  [scanned]" if (
            tuple(idxs) == q or tuple(idxs) == q[::-1]
        ) else ""
        lines.append(f"  {list(q)} ({types}): {'; '.join(pks)}{scanned}")
    return lines


def GetMultiDihedFcnFromIdxs(p,idxs):
    """ Get a MultiDihedFcn object from the Parm object p using the given indices.
    
    Parameters
    ----------
    p : parmed.AmberParm
        The Parm object containing the molecular structure.
    idxs : list of int
        A list of 4 integers representing the indices of the atoms forming the dihedral.
    
    Returns
    -------
    MultiDihedFcn
        A MultiDihedFcn object representing the dihedral function for the given indices.
    
    """
    xs = FindDihedrals(p,idxs)
    prims = []
    for x in xs:
        t = x.type
        prims.append( PrimDihedFcn(t.phi_k,t.phase,t.per) )
    return MultiDihedFcn(idxs,prims)


def ChangeParmFromMultiDihedFcn(p,fcn):
    """ Change parameters from a MultiDihedFcn object to a new Parm object. 
    
    Parameters
    ----------
    p : parmed.AmberParm
        The Parm object to be modified.
    fcn : MultiDihedFcn
        The MultiDihedFcn object containing the new dihedral parameters.
    
    Returns
    -------
    parmed.AmberParm
        A new Parm object with the dihedral parameters changed according to the MultiDihedFcn object.
        
    """
    out = CopyParm(p)
    scee = 1.2
    scnb = 2.0
    xs = []
    for prim in fcn.prim:
        per = prim.per
        ph = prim.phase
        fc = prim.fc
        x = DihedralType(fc, per, ph, scee, scnb)
        xs.append(x)
    ChangeDihedrals(out,fcn.idxs,xs)
    return out




# def GetDihedClasses(idxs=None):
    
#     class1 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,3)] ) ]
               
#     class2 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,  0,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,180,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,180,2)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1, 0,2)] ) ]

#     class3 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,  0,2),
#                           PrimDihedFcn(1,  0,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,180,2),
#                           PrimDihedFcn(1,180,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,180,2),
#                           PrimDihedFcn(1,180,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,  0,2),
#                           PrimDihedFcn(1,180,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,180,2),
#                           PrimDihedFcn(1,  0,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,  0,2),
#                           PrimDihedFcn(1,180,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
#                           PrimDihedFcn(1,180,2),
#                           PrimDihedFcn(1,  0,3)] ),
#               MultiDihedFcn( idxs, [PrimDihedFcn(1,180,1),
#                           PrimDihedFcn(1,  0,2),
#                           PrimDihedFcn(1,  0,3)] ) ]


#     return { 1: class1,
#              2: class2,
#              3: class3 }




def GetDihedClasses(idxs=None):
    """ Get a dictionary of MultiDihedFcn classes based on the provided indices.
    
    Parameters
    ----------
    idxs : list of int, optional
        A list of 4 integers representing the indices of the atoms forming the dihedral. If None, defaults to [0, 1, 2, 3].
    
    Returns
    -------
    dict
        A dictionary where keys are integers representing the number of terms in the dihedral function, and values are lists of MultiDihedFcn objects.
        
    """
    class1 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1)] ),
              MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,2)] ),
              MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,3)] ) ]


    class1 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1)] ) ]


    
    # class2 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
    #                                 PrimDihedFcn(1,  0,2)] ) ]

    # class3 = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,1),
    #                                 PrimDihedFcn(1,  0,2),
    #                                 PrimDihedFcn(1,  0,3)] ) ]

    cs = {}
    cs[1] = class1
    for j in range(2,13):
        cs[j] = [MultiDihedFcn( idxs, [PrimDihedFcn(1,  0,n) for n in range(1,j+1)])]
    

    return cs





def EnergyScansWithoutDihedrals(mol,list_of_los,cons):
    """ Calculate energies for a list of scans without dihedrals.
    
    Parameters
    ----------
    stdargs
        The standard arguments containing the molecular structure and calculator.
    list_of_scans : list of list of ase.Atoms
        A list of scans, where each scan is a list of ase.Atoms objects representing different geometries.
    cons : list of Constraint
        A list of Constraint objects representing the constraints to be applied to the geometries.
    
    Returns
    -------
    list of list of float
        A list of lists, where each inner list contains the energies for the corresponding scan.
    
    """
    from . AmberParm import CopyParm
    from . Dihedrals import DeleteDihedrals
    from . Dihedrals import GetMultiDihedFcnFromIdxs
    from tempfile import mkstemp
    import os,copy
    
    p = CopyParm(mol)
    llos = copy.deepcopy(list_of_los)
    
    DeleteDihedrals(p,[ x.idxs for x in cons ])

    fd,path = mkstemp(dir=".",prefix="tmp.",suffix=".parm7")
    if not os.isatty(fd):  # Check if fd is still valid
        os.close(fd)
    #fh = os.fdopen(fd,"w")
    #p.save("tmp.parm7",overwrite=True)
    p.save(path,overwrite=True,format="amber")
    for los in llos:
        los.calc = None
        for s in los:
            s.data["parm"] = path
    
    #calc = stdargs.MakeCalc(parm=path)
    
    list_of_enes = []
    for los in llos:
        enes = []
        
        for geom in los:
            t = geom.copy()
            t.constraints = None
            t.restraints = None
            t.data["constraints"] = []
            t.data["restraints"] = []
            g = t.GetASEAtoms()
            g.calc=los.BuildCalc(t)
            #g.calc = calc
            #g.calc.reset()
            e = g.get_potential_energy()
            enes.append(e)
        list_of_enes.append(enes)
        
    if os.path.exists(path):
        os.remove(path)
        
    for los in llos:
        los.calc = None
        
    return list_of_enes



def IsolatedLinearSolve(mol,idxs,losll,hlenes,nprim,pname, instance_idxs=None):
    """ Solve the isolated linear problem for dihedral parameters.
    
    Parameters
    ----------
    stdargs
        The standard arguments containing the molecular structure and calculator.
    idxs : list of int
        A list of 4 integers representing the indices of the atoms forming the dihedral.
    llgeoms : list of ase.Atoms
        A list of ase.Atoms objects representing the geometries for the low-energy scan.
    hlenes : list of float
        A list of floats representing the energies for the high-energy scan.
    nprim : int
        The number of primitive terms in the dihedral function.
    pname : str
        The name of the dihedral parameter set.
    
    Returns
    -------
    MultiDihedFcn
        A MultiDihedFcn object representing the best-fit dihedral function for the given parameters.
        
    """
    #from . AmberParm import GetDihedClasses
    from . Constraints import FillConstraints
    from . Constraints import Constraint
    import numpy as np
    import copy
    from . constants import AU_PER_KCAL_PER_MOL
    from . constants import AU_PER_ELECTRON_VOLT

    KCAL_PER_EV = AU_PER_ELECTRON_VOLT() / AU_PER_KCAL_PER_MOL()

    print(f"[fit] ---- IsolatedLinearSolve {pname} quartet={list(idxs)} nprim={nprim} ----")
    print(f"[fit] orig parm Fourier: {format_parm_quartet(mol, idxs)}")
    rotors = summarize_rotors_on_bond(mol, idxs)
    if rotors:
        print(f"[fit] proper terms on central bond {idxs[1]}-{idxs[2]}:")
        for line in rotors:
            print(f"[fit] {line}")

    inst = []
    seen = set()
    for q in list(instance_idxs or []) + [list(idxs)]:
        q = tuple(int(i) for i in q)
        key = q if q <= q[::-1] else q[::-1]
        if key in seen:
            continue
        seen.add(key)
        inst.append(list(q))
    if not inst:
        inst = [list(idxs)]
    print(
        f"[fit] leftover MM: delete {len(inst)} type-equivalent quartet(s) "
        f"{inst}  ({len(losll)} geometries, sander, no GeomOpt)"
    )

    graph = losll.structs[0].GetGraph()
    list_of_scans = [ losll ]
    scan_con = Constraint("dihed", idxs, graph=graph)
    cons = [ Constraint("dihed", q, graph=graph) for q in inst ]
    list_of_enes = EnergyScansWithoutDihedrals(mol,list_of_scans,cons)
    llenes = np.array(list_of_enes[0])
    hlenes = np.array(hlenes,copy=True)

    llenes *= KCAL_PER_EV
    hlenes *= KCAL_PER_EV

    # for igeom,llgeom in enumerate(losll):
    #     g = llgeom.GetASEAtoms()
    #     o = FillConstraints(g,cons,force=True)
    #     v = o[0].value
    #     #if abs(360-v) < 0.5:
    #     #    v = 0
    #     print("%3i %8.2f %20.6f %20.6f"%( igeom,v,llenes[igeom],hlenes[igeom] ))

    
    hlenes -= np.amin(hlenes)
    llenes -= np.amin(llenes)
    

    angs = []
    inst_rows = []
    for igeom,llgeom in enumerate(losll):
        g = llgeom.GetASEAtoms()
        o = FillConstraints(g,[scan_con],force=True)
        v = o[0].value
        if abs(360-v) < 0.01:
            v=0
        angs.append(v)
        row = []
        for q in inst:
            a = float(g.get_dihedral(*q))
            if abs(a - 360) < 0.01:
                a = 0.0
            row.append(a)
        inst_rows.append(row)

    data = []
    for i in range(len(losll)):
        data.append( [angs[i],llenes[i],hlenes[i], inst_rows[i]] )
    data = sorted(data,key=lambda x: x[0])
    angs = np.array( [x[0] for x in data] )
    llenes = np.array( [x[1] for x in data] )
    hlenes = np.array( [x[2] for x in data] )
    instance_angs = np.array( [x[3] for x in data], dtype=float )
    y = hlenes-llenes
    
    npts = len(y)
    y_ptp = float(np.max(y) - np.min(y)) if npts else 0.0
    y_std = float(np.std(y)) if npts else 0.0
    print(
        f"[fit] leftover y: n={npts} ptp={y_ptp:.3f} std={y_std:.3f} kcal/mol "
        f"φ=[{float(np.min(angs)):.1f}, {float(np.max(angs)):.1f}]°"
    )
    offsets = []
    for j in range(instance_angs.shape[1]):
        d = (instance_angs[:, j] - angs + 180.0) % 360.0 - 180.0
        offsets.append(float(np.mean(d)))
    print(
        f"[fit] {len(inst)} instance(s); mean φ offsets vs scan: "
        + ", ".join(f"{o:.1f}°" for o in offsets)
    )
    if len(inst) >= 3:
        print(
            f"[fit] {pname}: equivalent rotors on one bond — n=1,2 often cancel "
            "in the total energy; instance-sum LS should pick n=3 (or a multiple)"
        )

    from .DihedFitRegularize import (
        append_fit_trace,
        cap_effective_ptp,
        chemical_barrier_cap,
        clip_dihed_fcs,
        dense_torsion_ptp,
        effective_fourier,
        fit_leftover_fourier,
        format_prims,
        fourier_rss,
        phase_variant_functions,
        scale_fcs_to_ptp,
        _solve_leftover_linear,
    )
    import os as _os

    def _pad_to_nprim(dfcn):
        if dfcn is None:
            return MultiDihedFcn(
                idxs, [PrimDihedFcn(0.0, 0.0, n) for n in range(1, nprim + 1)]
            )
        if len(dfcn.prims) == nprim:
            return dfcn
        by_per = {}
        for p in dfcn.prims:
            key = (int(p.per), round(float(p.phase) % 360.0, 1))
            by_per[key] = by_per.get(key, 0.0) + float(p.fc)
        prims = []
        for n in range(1, nprim + 1):
            fc = 0.0
            phase = 0.0
            for (per, ph), val in by_per.items():
                if per == n:
                    fc = val
                    phase = ph
                    break
            prims.append(PrimDihedFcn(fc, phase, n))
        return MultiDihedFcn(idxs, prims)

    orig_dfcn = GetMultiDihedFcnFromIdxs(mol, idxs)
    rss_orig, _c_orig, _v_orig, r2_orig = fourier_rss(
        angs, y, orig_dfcn, instance_angs=instance_angs
    )
    print(
        f"[fit] orig vs leftover (instance-sum): {format_prims(orig_dfcn)}  "
        f"rss={rss_orig:.4g} r²={r2_orig:.4f}"
    )

    bestdfcn, _x, fit_info = fit_leftover_fourier(
        angs, y, nprim, idxs, pname=pname, instance_angs=instance_angs
    )
    phase_search = str(_os.environ.get("FFPOPT_DIHED_PHASE_SEARCH", "")).strip().lower() in {
        "1", "true", "yes", "on",
    }
    if phase_search:
        print(
            f"[fit] {pname}: extra 0°/180° phase search "
            "(FFPOPT_DIHED_PHASE_SEARCH=1; signed PK already covers this)"
        )
        best_rss, _, _, _ = fourier_rss(
            angs, y, bestdfcn, instance_angs=instance_angs
        )
        for dfcn in phase_variant_functions(idxs, len(bestdfcn.prims)):
            A = np.zeros((npts, len(dfcn.prims) + 1))
            for i, prim in enumerate(dfcn.prims):
                A[:, i] = np.sum(
                    np.cos(
                        np.deg2rad(
                            float(prim.per) * instance_angs + float(prim.phase)
                        )
                    ),
                    axis=1,
                )
            A[:, -1] = 1.0
            xph, _info = _solve_leftover_linear(A, y)
            dfcn.SetFCs(clip_dihed_fcs(xph[: len(dfcn.prims)]))
            rss_ph, _, _, _ = fourier_rss(
                angs, y, dfcn, instance_angs=instance_angs
            )
            if rss_ph < best_rss:
                best_rss = rss_ph
                bestdfcn = dfcn

    if bestdfcn is None:
        raise RuntimeError(f"IsolatedLinearSolve failed for {pname}")
    bestdfcn = _pad_to_nprim(bestdfcn)
    leftover_slack = 1.15
    try:
        leftover_slack = float(
            _os.environ.get("FFPOPT_DIHED_LEFTOVER_SLACK", "1.15")
        )
    except (TypeError, ValueError):
        leftover_slack = 1.15
    leftover_cap = max(y_ptp * leftover_slack, 0.5)
    chem_cap = chemical_barrier_cap(pname)
    target = min(leftover_cap, chem_cap)
    v_eff = effective_fourier(bestdfcn, instance_angs)
    ptp_eff = float(np.max(v_eff) - np.min(v_eff)) if v_eff.size else 0.0
    print(
        f"[fit] {pname}: one-quartet ptp={dense_torsion_ptp(bestdfcn):.2f}  "
        f"instance-sum ptp={ptp_eff:.2f} leftover ptp={y_ptp:.2f} "
        f"caps leftover={leftover_cap:.2f} chemical={chem_cap:.2f}"
    )
    cap_effective_ptp(bestdfcn, instance_angs, target, where=pname)
    q_ptp = dense_torsion_ptp(bestdfcn)
    if q_ptp > chem_cap:
        print(
            f"[fit] one-quartet ptp {q_ptp:.2f} > chemical cap {chem_cap:.2f} "
            f"at {pname} (cancelled harmonics leaked into a single oxygen). "
            "Scaling the written Fourier."
        )
        scale_fcs_to_ptp(bestdfcn, chem_cap)

    rss_fit, const_fit, v_fit, r2_fit = fourier_rss(
        angs, y, bestdfcn, instance_angs=instance_angs
    )
    min_rel = 0.02
    try:
        min_rel = float(_os.environ.get("FFPOPT_DIHED_MIN_REL_IMPROVE", "0.02"))
    except (TypeError, ValueError):
        min_rel = 0.02
    orig_ptp = dense_torsion_ptp(orig_dfcn) if orig_dfcn.prims else 0.0
    keep_orig = False
    if orig_ptp < 0.05:
        if r2_fit < 0.05:
            keep_orig = True
            print(
                f"[fit] KEEP ORIG at {pname}: scanned quartet is flat in GAFF "
                f"and leftover has no Fourier signal (r²={r2_fit:.4f}). "
                "orig vs itNN _dihed.png will stay flat/identical."
            )
    elif rss_fit >= rss_orig * (1.0 - min_rel):
        keep_orig = True
        print(
            f"[fit] KEEP ORIG at {pname}: fitted rss={rss_fit:.4g} "
            f"does not beat orig rss={rss_orig:.4g} "
            f"(need {100.0 * min_rel:.1f}% RSS drop). "
            "parm7 Fourier will stay GAFF — orig vs itNN _dihed.png will match."
        )
    if keep_orig:
        bestdfcn = _pad_to_nprim(copy.deepcopy(orig_dfcn))
        rss_fit, const_fit, v_fit, r2_fit = fourier_rss(
            angs, y, bestdfcn, instance_angs=instance_angs
        )

    v_eff = effective_fourier(bestdfcn, instance_angs)
    ptp_eff = float(np.max(v_eff) - np.min(v_eff)) if v_eff.size else 0.0
    print(
        f"[fit] {pname}: {format_prims(bestdfcn)}  "
        f"one-quartet ptp={dense_torsion_ptp(bestdfcn):.2f} "
        f"instance-sum ptp={ptp_eff:.2f} kcal/mol "
        f"rss={rss_fit:.4g} r²={r2_fit:.4f} keep_orig={keep_orig} "
        f"ninst={len(inst)}"
    )
    append_fit_trace(
        {
            "pname": pname,
            "idxs": list(idxs),
            "nprim": int(nprim),
            "y_ptp": y_ptp,
            "orig": format_prims(orig_dfcn),
            "fitted": format_prims(bestdfcn),
            "rss_orig": rss_orig,
            "rss_fit": rss_fit,
            "r2_orig": r2_orig,
            "r2_fit": r2_fit,
            "keep_orig": keep_orig,
            "lam": fit_info.get("lam"),
            "cond": fit_info.get("cond"),
            "nprim_sel": fit_info.get("nprim"),
            "n_irls": fit_info.get("n_irls"),
        }
    )

    fh = open(f"iso.{pname}.dat","w")
    fh.write("# %s\n"%(str(bestdfcn)))
    fh.write(f"# orig {format_prims(orig_dfcn)}\n")
    fh.write(f"# keep_orig={keep_orig} rss_orig={rss_orig:.6g} rss_fit={rss_fit:.4g}\n")
    for i in range(npts):
        fh.write("%12.3f %20.10e %20.10e %20.10e\n"%\
                 ( angs[i], hlenes[i], llenes[i],
                   llenes[i] + v_fit[i] + const_fit ) )
    fh.close()
    return bestdfcn
    


def AngularStdDev(angs):
    """ Calculate the circular standard deviation of a list of angles.
    
    Parameters
    ----------
    angs : list of float
        A list of angles in degrees.
    
    Returns
    -------
    float
        The circular standard deviation of the angles in degrees.
    
    """
    import numpy as np
    from scipy.stats import circstd
    rads = np.deg2rad(np.array(angs))
    return np.rad2deg(circstd(rads))



class ParamType(object):
    """ A class representing a type of dihedral parameter.
    
    Parameters
    ----------
    name : str
        The name of the dihedral parameter type.
    nprim : int
        The number of primitive terms in the dihedral function.
    masks : list of list of str
        A list of masks, where each mask is a list of 4 strings representing the atom
        indices forming the dihedral.
    
    Attributes
    ----------
    name : str
        The name of the dihedral parameter type.
    nprim : int
        The number of primitive terms in the dihedral function.
    masks : list of list of str
        A list of masks, where each mask is a list of 4 strings representing the atom
        indices forming the dihedral.
    dfcns : list of MultiDihedFcn, optional
        A list of MultiDihedFcn objects representing the dihedral functions for this parameter type.
    
    """
    
    def __init__(self,name,nprim,masks):
        self.name = name
        self.nprim = nprim
        self.masks = masks
        self.dfcns = None


        
    def get_dict(self):
        """ Get a dictionary representation of the ParamType object. """
        d = {}
        d[self.name] = { "nprim": self.nprim,
                         "masks": self.masks }
        return d

    
        
class ParamInstance(object):
    """ A class representing an instance of a dihedral parameter type.
    
    Parameters
    ----------
    ptype : ParamType
        The ParamType object representing the type of dihedral parameter.
    masks : list of list of str
        A list of masks, where each mask is a list of 4 strings representing the atom
        indices forming the dihedral.
    mol : parmed.AmberParm
        The Parm object containing the molecular structure.
    impropers : bool, optional
        If True, only consider improper dihedrals. Default is False.
    
    Attributes
    ----------
    ptype : ParamType
        The ParamType object representing the type of dihedral parameter.
    masks : list of list of str
        A list of masks, where each mask is a list of 4 strings representing the atom
        indices forming the dihedral.
    dihedidxs : list of tuple of int
        A list of tuples, where each tuple contains 4 integers representing the indices of the atoms
        forming the dihedral.

    """
    
    def __init__(self,ptype,masks,mol,impropers=False):
        from parmed.amber.mask import AmberMask
        
        self.ptype = ptype
        self.masks = masks
        self.dihedidxs = []
        
        possible_diheds = []
        for mask in masks:
            ats = []
            for i in range(4):
                sidxs = [idx for idx in AmberMask(mol,mask[i]).Selected()]
                ats.append(sidxs)
            for aidx in ats[0]:
                #amask = f"@{aidx+1}"
                for bidx in ats[1]:
                    #bmask = f"@{bidx+1}"
                    for cidx in ats[2]:
                        #cmask = f"@{cidx+1}"
                        for didx in ats[3]:
                            #dmask = f"@{didx+1}"
                            fwd = (aidx,bidx,cidx,didx)
                            possible_diheds.append(fwd)
            for tgt in possible_diheds:
                for d in mol.dihedrals:
                    if d.improper and not impropers:
                        continue
                    fwd = (d.atom1.idx,d.atom2.idx,d.atom3.idx,d.atom4.idx)
                    rev = tuple(list(fwd)[::-1])
                    if fwd == tgt or rev == tgt:
                        if tgt not in self.dihedidxs:
                            self.dihedidxs.append(tgt)
                            


    def get_dict(self):
        """ Get a dictionary representation of the ParamInstance object. """
        from collections import defaultdict as ddict
        d = ddict(list)
        name = self.ptype.name
        for didx in self.dihedidxs:
            d[name].append( [ f"@{didx[0]+1}",f"@{didx[1]+1}",
                              f"@{didx[2]+1}",f"@{didx[3]+1}" ] )
        return d


    
class ProfileType(object):
    """ A class representing a profile for dihedral parameters.
    
    Parameters
    ----------
    hl : str
        The name of the high-energy scan file.
    ll : str
        The name of the low-energy scan file.
    name : str
        The name of the profile.
    plots : list of str
        A list of strings representing the plots to be generated for this profile.
    stdargs
        The standard arguments containing the molecular structure and calculator.
    
    Attributes
    ----------
    hl : str
        The name of the high-energy scan file.
    ll : str
        The name of the low-energy scan file.
    name : str
        The name of the profile.
    plots : list of str
        A list of strings representing the plots to be generated for this profile.
    #hlatoms : list of ase.Atoms
    #    A list of ase.Atoms objects representing the high-energy scan geometries.
    #hlcons : list of Constraint
    #    A list of Constraint objects representing the constraints for the high-energy scan geometries.
    #llatoms : list of ase.Atoms
    #    A list of ase.Atoms objects representing the low-energy scan geometries.
    #llcons : list of Constraint
    #    A list of Constraint objects representing the constraints for the low-energy scan geometries.
    loshl : ffpopt.Struct.ListOfStruct
         The list of high-level calculations
    losll : ffpopt.Struct.ListOfStruct
         The list of low-level calculations
    """
    def __init__(self,hl,ll,name,plots,stride):
        #from . Reader import ReadGeomsFromXYZ
        from . Struct import ListOfStruct
        self.hl = hl
        self.ll = ll
        self.name = name
        self.plots = plots
        #self.hlatoms,self.hlcons = ReadGeomsFromXYZ(stdargs,self.hl)
        #self.llatoms,self.llcons = ReadGeomsFromXYZ(stdargs,self.ll)
        self.loshl = ListOfStruct.from_file(self.hl)
        self.losll = ListOfStruct.from_file(self.ll)

        
        if len(self.loshl.structs) != len(self.losll.structs):
            raise Exception(f"Structure count mismatch in {self.hl} "
                            +f"and {self.ll} ({len(self.loshl)} vs "
                            +f"{len(self.losll)})")
        
        s = stride
        self.loshl.structs = self.loshl.structs[::s]
        self.losll.structs = self.losll.structs[::s]

        #self.hlatoms = self.hlatoms[::s]
        #self.hlcons = self.hlcons[::s]
        #self.llatoms = self.llatoms[::s]
        #self.llcons = self.llcons[::s]

        
    def get_dict(self):
        """ Get a dictionary representation of the ProfileType object. """
        return { "hl": self.hl, "ll": self.ll,
                 "name": self.name, "plots": self.plots }

        
        
class SystemType(object):
    """ A class representing a system for dihedral parameter fitting.
    
    Parameters
    ----------
    parm : Amber parm7 filename
        Original parameter file
    output : str
        The name of the output file for the fitted parameters.
    profilelist : list of dict
        A list of dictionaries representing the profiles for this system.
    pdict : dict
        A dictionary where keys are parameter names and values are lists of masks for those parameters.
    ptypedict : dict
        A dictionary where keys are parameter names and values are ParamType objects representing the types of di
        hedral parameters.
    
    Attributes
    ----------
    parm : str, Amber parm7 filename
        The filename of the original parm7 file
    mol : parmed parm7 structure
        The parmed representation of the parm7 file
    output : str
        The name of the output file for the fitted parameters.
    pinstances : list of ParamInstance
        A list of ParamInstance objects representing the instances of dihedral parameters for this system.
    profiles : list of ProfileType
        A list of ProfileType objects representing the profiles for this system.
    
    """
    
    def __init__(self,parm,output,profilelist,pdict,ptypedict,stride=1):
        import parmed
        #self.stdargs = stdargs
        self.parm = parm
        self.mol = parmed.load_file(self.parm)
        self.output = output
        self.pinstances = []
        self.profiles = []
        
        for pname in pdict:
            if pname not in ptypedict:    
                raise Exception(f"Parameter {pname} in {parm} "
                                +f"not found in parameter list "
                                +f"{[name for name in ptypedict]}")
            
            ptype = ptypedict[pname]
            masks = pdict[pname]
            self.pinstances.append( ParamInstance(ptype,masks,self.mol) )

        for pname in ptypedict:
            if pname not in pdict:
                ptype = ptypedict[pname]
                masks = ptype.masks
                pinst = ParamInstance(ptype,masks,self.mol)
                if len(pinst.dihedidxs) > 0:
                    self.pinstances.append(pinst)
            
        for profile in profilelist:
            for key in ["hl","ll","name","plots"]:
                if key not in profile:
                    raise Exception(f"'profiles' section is missing '{key}' key")
            hl = profile["hl"]
            ll = profile["ll"]
            name = profile["name"]
            plots = profile["plots"]
            self.profiles.append( ProfileType(hl,ll,name,plots,stride) )

            

    def make_new_parm(self):
        """ Create a new Parm object with the dihedral parameters set according to the instances.
        
        Returns
        -------
        parmed.AmberParm
            A new Parm object with the dihedral parameters set according to the instances.
        
        """
        from parmed import DihedralType
        from . Dihedrals import ChangeDihedrals
        from . AmberParm import CopyParm

        p = CopyParm(self.mol)
        
        scee = 1.2
        scnb = 2.0
        for pinst in self.pinstances:
            xs = []
            for prim in pinst.ptype.dfcns.prims:
                per = prim.per
                ph = prim.phase
                fc = prim.fc
                x = DihedralType(fc, per, ph, scee, scnb)
                xs.append(x)
            for idxs in pinst.dihedidxs:
                ChangeDihedrals(p,idxs,xs)
        return p

    
            
    def get_dict(self):
        """ Get a dictionary representation of the SystemType object.
        
        Returns
        -------
        dict
            A dictionary containing the standard arguments, output file name, parameter instances, and profiles.
        
        """
        params = {}
        profiles = []

        for p in self.pinstances:
            x = p.get_dict()
            for key in x:
                params[key] = x[key]

        for p in self.profiles:
            profiles.append(p.get_dict())
        
        d = { "parm": self.parm,
              #"crd": self.stdargs.crd,
              "output": self.output,
              "params": params,
              "profiles": profiles }
        return d


    
    def find_pinstance(self,pname):
        """ Find a ParamInstance by its name.
        
        Parameters
        ----------
        pname : str
            The name of the parameter type to find.
        
        Returns
        -------
        ParamInstance or None
            The ParamInstance object if found, otherwise None.
            
        """
        p=None
        for pinst in self.pinstances:
            if pinst.ptype.name == pname:
                p = pinst
        return p


    def write_output(self):
        """ Write the output file with the dihedral parameters.
        
        This method creates a new Parm object with the dihedral parameters set according to the instances,
        and writes the parameters to the specified output file using the WriteParmedScript function.
        
        Returns
        -------
        None
        
        """
        import copy
        dfcns = []
        for pinst in self.pinstances:
            for idxs in pinst.dihedidxs:
                dfcn = copy.deepcopy(pinst.ptype.dfcns)
                dfcn.idxs = idxs
                dfcns.append(dfcn)
        WriteParmedScript(self.output,self.mol,dfcns)
    
    
class FitInputType(object):
    """ A class representing the input for fitting dihedral parameters.
    
    This class is used to read the input data from a JSON file and create instances of ParamType, SystemType, and ProfileType.
    It also provides methods to create initial guesses for the dihedral parameters and to write the output in a specific format.
    
    Parameters
    ----------
    args
        The standard arguments containing the molecular structure and calculator.
    datadict : dict
        A dictionary containing the input data for fitting dihedral parameters.
    
    Attributes
    ----------
    iteration : int
        The current iteration number for fitting dihedral parameters.
    ptypedict : dict
        A dictionary where keys are parameter names and values are ParamType objects representing the types of di
        hedral parameters.
    output : str
        The name of the output file for the fitted parameters.
    systems : list of SystemType
        A list of SystemType objects representing the systems for which dihedral parameters are to be fitted
    
    """
    
    @classmethod
    def from_file(cls,args,fname):
        import json
        data = {}
        with open(fname,"r") as file:
            data = json.load(file)
        return cls(args,data)
    
        
    def __init__(self,args,datadict):
        #from . Options import StandardArgs
        self.iteration = 0

        for key in ["params","output","systems"]:
            if key not in datadict:
                raise Exception(f"Input is missing '{key}' key")
                
        self.ptypedict = {}
        for pname in datadict["params"]:
            nprim = datadict["params"][pname]["nprim"]
            masks = datadict["params"][pname]["masks"]
            if masks is not None:
                for mask in masks:
                    for atmask in mask:
                        if atmask[:2] != "@%":
                            raise Exception(f"Global dihedral {mask} contains mask {atmask}, which should start with @%%")
            self.ptypedict[pname] = ParamType(pname,nprim,masks)
            #print(pname)
        self.output = datadict["output"]
        self.systems = []
        for s in datadict["systems"]:
            for key in ["parm","output","params","profiles"]:
                if key not in s:
                    raise Exception(f"'systems' section is missing '{key}' key")
            parm = s["parm"]
            #crd = s["crd"]
            output = s["output"]
            pdict = s["params"]
            profilelist = s["profiles"]
            #stdargs = StandardArgs.from_manual_parm(parm,crd,args)
            self.systems.append\
                ( SystemType\
                  (parm,output,profilelist,pdict,self.ptypedict) )

        unused = []
        for pname in self.ptypedict:
            found=False
            for s in self.systems:
                for pinst in s.pinstances:
                    if pinst.ptype.name == pname:
                        found=True
            if not found:
                unused.append(pname)
        for pname in unused:
            print(f"Removing parameter {pname} because it was unused")
            del self.ptypedict[pname]
                
            
    def get_dict(self):
        """ Get a dictionary representation of the FitInputType object.
        
        A dictionary containing the parameters, output file name, and systems.

        Returns
        -------
        dict
            A dictionary containing the parameters, output file name, and systems.

        """
        d = {}
        d["params"] = {}
        for pname in self.ptypedict:
            x = self.ptypedict[pname].get_dict()
            d["params"][pname] = x[pname]
        d["output"] = self.output
        d["systems"] = []
        for s in self.systems:
            d["systems"].append( s.get_dict() )

        return d

    def get_json(self):
        """ Get a JSON representation of the FitInputType object."""
        import json
        d = self.get_dict()
        return json.dumps(d,indent=4)
            

    def make_initial_guesses(self):
        """ Create initial guesses for the dihedral parameters based on the profiles.

        This method iterates through the parameter types and their instances, finds the best profile for each dihedral type,
        and computes the isolated linear solve to determine the dihedral function coefficients.

        Returns
        -------
        list
            A list of dihedral parameters computed from the best profiles for each parameter type.

        """

        for pname in self.ptypedict:
            #print(pname)
            bests = None
            bestpinst = None
            bestprof = None
            bestidxs = None
            beststd = -1
            
            for s in self.systems:
                for pinst in s.pinstances:
                    if pinst.ptype.name == pname:
                        for idxs in pinst.dihedidxs:
                            for prof in s.profiles:
                                angs = []
                                for struct in prof.losll:
                                    atoms = struct.GetASEAtoms()
                                    ang = atoms.get_dihedral(*idxs)
                                    if abs(ang-360) < 0.01:
                                        ang = 0.
                                    angs.append(ang)
                                    #print("ang=%.2f"%(ang))
                                if len(angs) > 2:
                                    astd = AngularStdDev(angs)
                                    #print("astd=%.2f %.2f"%(astd,beststd))
                                    if astd > beststd:
                                        beststd = astd
                                        bestprof = prof
                                        bestpinst = pinst
                                        bestidxs = idxs
                                        bests = s
            llgeoms = bestprof.losll
            hlenes = [ hlgeom.data["energy"]
                       for hlgeom in bestprof.loshl ]
            nprim = bestpinst.ptype.nprim

            
            dfcns = IsolatedLinearSolve\
                (bests.mol,bestidxs,llgeoms,hlenes,nprim,pname,
                 instance_idxs=bestpinst.dihedidxs)
            
            bestpinst.ptype.dfcns = dfcns
        return self.get_params()

            
    def get_num_params(self):
        """ Get the total number of primitive parameters across all parameter types.
        
        Returns
        -------
        int
            The total number of primitive parameters across all parameter types.

        """

        n = 0
        for pname in self.ptypedict:
            n += self.ptypedict[pname].nprim
        return n

    

    def get_params(self):
        """ Get the current values of the dihedral function coefficients.
        
        Returns
        -------
        numpy.ndarray
            A numpy array containing the current values of the dihedral function coefficients for all parameter types.
        
        """
        import numpy as np
        x=[]
        for pname in self.ptypedict:
            ptype = self.ptypedict[pname]
            dfcns = ptype.dfcns
            for prim in dfcns.prims:
                x.append( prim.fc )
        return np.array(x)


    
    def set_params(self,x):
        """ Set the dihedral function coefficients based on the provided array.
        
        Parameters
        ----------
        x : numpy.ndarray
            A numpy array containing the new values for the dihedral function coefficients.
        
        """
        ipar = 0
        for pname in self.ptypedict:
            ptype = self.ptypedict[pname]
            nprim = ptype.nprim
            ptype.dfcns.SetFCs( x[ipar:ipar+nprim] )
            ipar += nprim

    def write_output(self):
        """ Write the output file with the dihedral parameters.
        
        This method creates a new AmberParameterSet object, populates it with the dihedral types from the parameter types,
        and writes the parameters to the specified output file.
        
        Returns
        -------
        None
        
        """
        from parmed.amber import AmberParameterSet
        from parmed import DihedralTypeList
        from parmed import DihedralType
        #from parmed.amber.mask import AmberMask
        
        for s in self.systems:
            s.write_output()
            
        fmod = AmberParameterSet()
        for pname in self.ptypedict:
            ptype = self.ptypedict[pname]
            if ptype.masks is not None and ptype.dfcns is not None:

                typs = DihedralTypeList()
                scee = 1.2
                scnb = 2.0
                for iprim,prim in enumerate(ptype.dfcns.prims):
                    per = prim.per
                    ph = prim.phase
                    fc = prim.fc
                    typ = DihedralType(fc, per, ph, scee, scnb)
                    typs.append(typ)
                
                for mask in ptype.masks:
                    atypes = [ x[2:] for x in mask ]
                    fwd = tuple(atypes)
                    rev = tuple(list(fwd)[::-1])
                    fmod.dihedral_types[fwd] = typs
                    fmod.dihedral_types[rev] = typs
                    
        fmod.write(self.output)
    

def DihedFitObjFcn(x,self):
    """ Objective function for fitting dihedral parameters.
    
    This function computes the chi-squared value based on the dihedral parameters
    and the energies of the low-energy and high-energy scans. It updates the positions
    of the low-energy geometries based on the computed energies and writes the results
    to files for each profile in the system.
    
    Parameters
    ----------
    x : numpy.ndarray
        A numpy array containing the dihedral function coefficients to be optimized.
    self : FitInputType
        An instance of FitInputType containing the systems and profiles for fitting.
    Returns
    ------- 
    float
        The chi-squared value representing the difference between the computed and expected energies.
    
    """
    import numpy as np
    from . GeomOpt import GeomOpt
    from . constants import AU_PER_KCAL_PER_MOL
    from . constants import AU_PER_ELECTRON_VOLT
    from tempfile import mkstemp
    import os
    
    KCAL_PER_EV = AU_PER_ELECTRON_VOLT() / AU_PER_KCAL_PER_MOL()
    from .DihedFitRegularize import dihed_fc_abs_max

    cap = dihed_fc_abs_max()
    x = np.clip(np.asarray(x, dtype=float), -cap, cap)
  
    chisq = 0

    it = self.iteration
    
    self.set_params(x)
    for isys,s in enumerate(self.systems):
        p = s.make_new_parm()
        
        fd,path = mkstemp(dir=".",prefix="tmp.",suffix=".parm7")
        if not os.isatty(fd):  # Check if fd is still valid
            os.close(fd)
        
        #p.save("tmp.parm7",overwrite=True)
        #s.stdargs.calc = s.stdargs.MakeCalc(parm="tmp.parm7")

        p.save(path,overwrite=True,format="amber")    
        #s.stdargs.calc = s.stdargs.MakeCalc(parm=path)
        #if os.path.exists(path):
        #    os.remove(path)

        

        

        for iprof,prof in enumerate(s.profiles):
            llene = []
            hlene = []
            prof.losll.calc = None

            for igeom in range(len(prof.losll)):
                #atoms = prof.llatoms[igeom].copy()
                #cons = prof.llcons[igeom]
                #print("cons=",igeom,cons)
                #atoms = GeomOpt(atoms,s.stdargs,constraints=cons)
                inpstruct = prof.losll.structs[igeom].copy()
                inpstruct.data["parm"] = path
                sout = GeomOpt( prof.losll, inpstruct )
                hlene.append( prof.loshl.structs[igeom].data["energy"] * KCAL_PER_EV )
                llene.append( sout.data["energy"] * KCAL_PER_EV )
                prof.losll.structs[igeom].Update\
                    ( sout.get_potential_energy(),
                      sout.get_positions(),
                      sout.get_forces())
                #prof.llatoms[igeom].set_positions( atoms.get_positions() )
                
            llene  = np.array(llene)
            llene -= np.amin(llene)
            
            hlene  = np.array(hlene)
            hlene -= np.amin(hlene)
            
            d = hlene-llene
            llene += np.mean(d)
            d = hlene-llene
            mychisq = np.dot(d,d)
            
            chisq += mychisq
            
            hlene -= np.amin(hlene)
            llene -= np.amin(llene)

            # Per-COBYLA-step dumps (mfit.itXX.*.NNNN.dat). Off by default:
            # nlmaxiter=300 × bonds × fragments floods fragment libraries.
            # Re-enable with FFPOPT_WRITE_MFIT=1 for ffpopt-DihedTwistAnimate --source mfit.
            write_mfit = os.environ.get("FFPOPT_WRITE_MFIT", "").strip().lower() in {
                "1", "true", "yes", "on",
            }
            if not write_mfit:
                continue
            if prof.name is None or prof.plots is None:
                continue
            elif len(prof.plots) == 0:
                continue

            for pname in prof.plots:
                pinst = s.find_pinstance(pname)
                if pinst is None:
                    continue
                            
                for idxs in pinst.dihedidxs:
                    angs = []
                    for struct in prof.losll:
                        atoms = struct.GetASEAtoms()
                        ang = atoms.get_dihedral( *idxs )
                        angs.append(ang)
                    data = []
                    for i in range(len(angs)):
                        data.append( [angs[i],hlene[i],llene[i]] )
                        #print("ang=",i,angs[i])
                    data = sorted(data,key=lambda x: x[0])
                
                    idxsname = "-".join( [ f"{i}" for i in idxs] )
                    fname = f"mfit.{prof.name}.{idxsname}.{it:04d}.dat"
                    fh = open(fname,"w")
                    fh.write("# %25.14f\n"%(mychisq))
                    for x in data:
                        fh.write("%20.10e %20.10e %20.10e\n"%(x[0],x[1],x[2]))
                    fh.close()
        if os.path.exists(path):
            os.remove(path)

    self.iteration += 1
    return chisq

                


def NonlinearSolve(args,finp):
    """ Perform a nonlinear optimization to fit dihedral parameters.
    
    Parameters
    ----------
    args : argparse.Namespace
        The command-line arguments containing the optimization parameters.
    finp : FitInputType
        An instance of FitInputType containing the systems and profiles for fitting.

    Returns
    -------
    None
        The function modifies the FitInputType instance in place, setting the optimized dihedral parameters.
    
    """
    from scipy.optimize import minimize
    import numpy as np
    import copy

    #objfcn = NonlinearObjective(stdargs,udscans)
    
    from .DihedFitRegularize import clip_dihed_fcs, dihed_fc_abs_max

    n = finp.get_num_params()
    x = finp.make_initial_guesses()
    cap = dihed_fc_abs_max()
    x = clip_dihed_fcs(x, where="linear-x0")
    finp.set_params(x)

    # Inner COBYLA+GeomOpt matches *relaxed total* MM energy and can walk
    # PKs up to the clip even when leftover is small (11 → 25 → 35 kcal).
    # The outer twist loop already rescans; default is the isolated linear fit.
    import os as _os

    relax = str(_os.environ.get("FFPOPT_FIT_RELAX", "")).strip().lower() in {
        "1", "true", "yes", "on",
    }
    if not relax:
        print(
            "[fit] regularized linear FCs only "
            "(set FFPOPT_FIT_RELAX=1 for COBYLA + per-step GeomOpt)"
        )
        return

    bounds = [(-cap, cap) for _ in range(n)]

    for s in finp.systems:
        for p in s.profiles:
            p.losll.SetArgs(args)

    res = minimize(
        DihedFitObjFcn,
        x,
        args=(finp,),
        method="COBYLA",
        bounds=bounds,
        options={
            "rhobeg": args.nlrhobeg,
            "tol": args.nltol,
            "maxiter": args.nlmaxiter,
            "disp": True,
        },
    )

    print(res)
    finp.set_params(clip_dihed_fcs(res.x, where="nonlinear-x"))



def WriteParmedScript(fname,p,dfcns): #,bytype):
    """ Write a Parmed script to modify dihedral parameters in a Parm object.
    
    Parameters
    ----------
    fname : str
        The name of the output file where the script will be written.
    p : parmed.AmberParm
        The Parm object containing the molecular structure.
    dfcns : list of MultiDihedFcn
        A list of MultiDihedFcn objects representing the dihedral functions to be modified.
    
    Returns
    -------
    None
    
    """
    from collections import defaultdict as ddict
    from . Dihedrals import FindDihedrals
    
    aidxs = [ idx for dfcn in dfcns for idx in dfcn.idxs ]
    aidxs = list(set(aidxs))
    resname = p.atoms[aidxs[0]].residue.name
    
    fh = open(fname,"w")
    fh.write("#!/usr/bin/env python3\n")
    fh.write("import sys\n")
    fh.write("import argparse\n")
    fh.write("from parmed import load_file\n")
    fh.write("from parmed.tools.actions import deleteDihedral, addDihedral\n")
    fh.write("from parmed.amber.mask import AmberMask\n")

    fh.write("parser = argparse.ArgumentParser(\"replace dihedral parameters\")\n")
    fh.write(f"parser.add_argument(\"--resname\",default=\"{resname}\",help=\"Name of the residue, default: {resname}\",type=str)\n")
    fh.write("parser.add_argument(\"iparm\",help=\"Input parm7\")\n")
    fh.write("parser.add_argument(\"oparm\",help=\"Output parm7\")\n")
    fh.write("args = parser.parse_args()\n")
    fh.write("rname = args.resname\n")

    fh.write("scee = 1.2\n")
    fh.write("scnb = 2.0\n")

    
    fh.write("if args.iparm == args.oparm:\n")
    fh.write("    raise Exception(\"The 2 filenames must be different\")\n\n")
    
    fh.write("p = load_file( args.iparm )\n\n")
    

    #if not bytype:
    if True:
        for aidx in aidxs:
            res = p.atoms[aidx].residue.name
            name = p.atoms[aidx].name
            mask = ":%s@%s"%("{rname}",name)
            fh.write(f"mask=f\"{mask}\"\n")
            fh.write("res = [i for i in AmberMask(p,mask).Selected()]\n")
            fh.write("if len(res) == 0:\n")
            fh.write("    raise Exception(f\"No atoms matching {mask}\")\n")
            

    n_replace = 0
    n_skip = 0
    fh.write("\n\n")
    for dfcn in dfcns:
        allmasks = [ [ ":%s@%s"%("{rname}",p.atoms[idx].name)
                       for idx in dfcn.idxs ] ]

        for masks in allmasks:
            mstr = ",".join(["f\"%s\""%(mask) for mask in masks])
            nonzero = [prim for prim in dfcn.prims if abs(float(prim.fc)) > 1.0e-8]
            orig_line = format_parm_quartet(p, dfcn.idxs)
            if not nonzero:
                n_skip += 1
                fh.write(
                    f"# skip quartet {dfcn.idxs}: all PKs ~ 0 "
                    "(keep original GAFF terms)\n\n"
                )
                print(
                    f"[fit] script {fname}: SKIP {list(dfcn.idxs)} "
                    f"(fitted PKs all ~0; orig {orig_line})"
                )
                continue
            n_replace += 1
            pk_str = ", ".join(
                f"n={int(prim.per)} PK={float(prim.fc):.4f} γ={float(prim.phase):.0f}"
                for prim in nonzero
            )
            print(
                f"[fit] script {fname}: REPLACE {list(dfcn.idxs)} "
                f"{pk_str}  (was {orig_line})"
            )
            fh.write(f"deleteDihedral(p,{mstr}).execute()\n")
            for prim in nonzero:
                fh.write(
                    f"addDihedral(p,{mstr},{prim.fc},{prim.per},{prim.phase},scee,scnb).execute()\n"
                )
            fh.write("\n\n")

    fh.write("p.save(args.oparm,overwrite=True)\n")
    fh.close()
    print(
        f"[fit] WriteParmedScript {fname}: replace={n_replace} skip_zero={n_skip}"
    )
    if n_replace == 0:
        print(
            f"[fit] WARNING {fname}: no quartet was replaced. "
            "itNN.parm7 Fourier == origparm; orig vs itNN _dihed.png will match."
        )

    







    
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################




def FindPuckers(s):
    """ Guess the C1,C2,C3,C4,O4 sugar atoms for each 5-membered ring in
    the structure

    Parameters
    ----------
    s : ffpopt.Struct.Struct
        The input structure to examine

    Returns
    -------
    rings : list of tuple
        The length of rings is the number of 5-membered rings
        Each element of the list is a tuple.
        The first element of the tuple is a list: the 5 indexes of the ring
        The second element is the same list sorted in order [C1,C2,C3,C4,O4]
        If the order could not be guessed, then the second element is None
    """
    
    g = s.GetGraph()
    mincycs = g.FindMinCycles()
    keepcycs = []
    for c in mincycs:
        if len(c) == 6:
            h = [ x for x in c[:-1] ]
            keepcycs.append(h)

    rings = []
    for c in keepcycs:
        ringseq = PuckerGuessByName(c,g,s)
        if ringseq is None:
            ringseq = PuckerGuessByElement(c,g,s)
        rings.append( (c,ringseq) )
        
    return rings


def PuckerGuessByName(cinp,g,s):
    """ Guess the C1,C2,C3,C4,O4 sugar atoms by element

    1. It assumes cinp is a list of 5 integers.
    2. The atom names are checked.
       2a. One atom name must contain "1"
       2b. One atom name must contain "2"
       2c. One atom name must contain "3"
       2d. Two atoms must contain "4"
    3. The O4 position is chosen from 2d by checking for
       a covalent bond to 2a.
    """
    from collections import defaultdict as ddict
    c = [ int(x) for x in cinp ]
    onames = [x for x in s.data["names"]]

    bad = False
    
    ignores = [ str(x) for x in range(10,100) ]
    for i in range(len(onames)):
        for ig in ignores:
            if ig in onames[i]:
                if i in c:
                    # One of the atoms in the ring has a name like
                    # C21 rather than C1 or C2 so we really can't
                    # use the names to figure out what is going on.
                    bad = True
                    
    dpos = ddict(list)
    unknown = []
    if not bad:
        #n = len(onames)
        #cmask = [False]*n
        #for i in c:
        #    cmask[i] = True
        for i in c:
            found=False
            for ipos in [1,2,3,4]:
                pos=str(ipos)
                if pos in onames[i]:
                    dpos[ ipos ].append( i )
                    found=True
                    break
            if not found:
                unknown.append(i)
            
    bad = len(unknown) > 0
    if not bad:
        if len(dpos[4]) != 2:
            bad=True
        elif len(dpos[3]) != 1:
            bad=True
        elif len(dpos[2]) != 1:
            bad=True
        elif len(dpos[1]) != 1:
            bad=True
            
    #hpos = ddict(list)
    #for x in dpos:
    #    for u in dpos[x]:
    #        hpos[x].append( onames[u] )

    if not bad:
        if str(dpos[1][0]) in g.edges[ str(dpos[4][0]) ]:
            O4 = dpos[4][0]
            C4 = dpos[4][1]
        elif str(dpos[1][0]) in g.edges[ str(dpos[4][1]) ]:
            O4 = dpos[4][1]
            C4 = dpos[4][0]
        else:
            bad = True

    ringseq = None
    if not bad:
        C1 = dpos[1][0]
        C2 = dpos[2][0]
        C3 = dpos[3][0]
        ringseq = [ C1,C2,C3,C4,O4 ]
    return ringseq




def PuckerGuessByElement(cinp,g,s):
    """
    Guess the C1,C2,C3,C4,O4 sugar atoms by element

    1. It assumes cinp is a list of 5 integers.
    2. 4 of the positions must correspond to "C", the non-carbon
       takes the O4' position.
    3. The bonding pattern is used to define a "clock". The
       clock either corresponds to O4,C1,C2,C3,C4 or
       O4,C4,C3,C2,C1 -- the rest of the algorithm is to figure
       out which of these is correct
    4. The C2/C3 positions are checked for bonded oxygens.
       a. Algorithm fails if either has more than 1 oxygen.
       b. It one has a O and the other doesn't, then it is assumed to be DNA, and the C3 position is chosen to be the one with the oxygen
       c. If they each have 1 oxygen, then the O4-C1-C2-O2 and O4-C4-C3-O3 dihedrals are calculated, and the decision is made based on these dihedrals.
    
    """

    
    from collections import defaultdict as ddict
    c = [ int(x) for x in cinp ]
    onames = [x for x in s.data["elements"]]

    nonc = [x for x in c if onames[x] != "C"]
    bad = len(nonc) != 1

    bonds = ddict(list)
    for a in g.edges:
        for b in g.edges[a]:
            bonds[ int(a) ].append( int(b) )

    clock = []
    if not bad:
        O4 = nonc[0]
        n1or4 = []
        for x in c:
            if x in bonds[O4]:
                n1or4.append(x)
        bad = len(n1or4) != 2
        
    if not bad:
        excl = n1or4 + [O4]
        clock = [O4,n1or4[0],None,None,n1or4[1]]
        for x in c:
            if x in excl:
                continue
            elif x in bonds[n1or4[0]]:
                clock[2] = x
            elif x in bonds[n1or4[1]]:
                clock[3] = x
        bad = None in clock
        
    if not bad:
        nbors2 = [ x for x in bonds[clock[2]] if x not in clock ]
        nnams2 = [ onames[x] for x in nbors2 ]
        nisO2 = [ 1 if name == "O" else 0 for name in nnams2 ]
        
        nbors3 = [ x for x in bonds[clock[3]] if x not in clock ]
        nnams3 = [ onames[x] for x in nbors3 ]
        nisO3 = [ 1 if name == "O" else 0 for name in nnams3 ]
        
        if sum(nisO2) > 1:
            bad = True
        elif sum(nisO3) > 1:
            bad = True

        if bad:
            clock=None
            return clock
        
        Oin2 = "O" in nnams2
        Oin3 = "O" in nnams3
        if Oin2 and not Oin3:
            # DNA
            clock = clock[::-1]
        elif Oin3 and not Oin2:
            # DNA
            pass
        elif not Oin2 and not Oin3:
            # Unidentifiable from this algorithm
            clock = None
            bad = True
        else:
            # RNA
            O2 = nbors2[ nnams2.index("O") ]
            #print(clock[0],clock[1],clock[2], O2)
            a2 = s.get_dihedral( clock[0],clock[1],clock[2], O2)
            #print(a2)
            O3 = nbors3[ nnams3.index("O") ]
            #print(clock[0],clock[4],clock[3], O3)
            a3 = s.get_dihedral( clock[0],clock[4],clock[3], O3)
            #print(a3)
            if a2 < 180 and a3 > 180:
                O4 = clock[0]
                C1 = clock[1]
                C2 = clock[2]
                C3 = clock[3]
                C4 = clock[4]
                clock = [C1,C2,C3,C4,O4]
                #pass
            elif a2 > 180 and a3 < 180:
                O4 = clock[0]
                C1 = clock[4]
                C2 = clock[3]
                C3 = clock[2]
                C4 = clock[1]
                clock = [C1,C2,C3,C4,O4]
                #clock = clock[::-1]
                pass
            else:
                bad = True
                clock = None
          
    return clock



