#!/usr/bin/env python3


def GeomOpt_ASE(los,struct,constraints=None,restraints=None):
    """
    Perform a structural geometry relaxation using the ASE BFGS minimization algorithm.

    This function synchronizes the input geometry to an ASE Atoms object, establishes
    any user-defined constraints or restraints by merging active lists with the original
    molecular parameters, binds the composite potential evaluator, and relaxes the atomic
    positions until convergence tolerances are satisfied.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The orchestrating collection tracking active command-line parameters and arguments.
    struct : ffpopt.Struct.Struct
        The initial input molecular geometry and parameter context to optimize.
    constraints : list of Constraint, optional
        A sequence of geometric constraints to merge and apply during minimization.
    restraints : list of Restraint, optional
        A sequence of structural restraints to merge and apply during minimization.

    Returns
    -------
    ffpopt.Struct.Struct
        The relaxed molecular structure with updated coordinates, energy, and forces.
    """

    from . CpuThreads import pin_math_threads
    pin_math_threads(1)
 
    import copy
    import ase.io
    from . Constraints import ConstraintList
    from . Restraints import RestraintList
    from . Constraints import ApplyConstraints, to_ase
    from . AseEngine import get_persistent_calc
    
    if True:
        
        from ase.optimize import BFGS

        reslist = None
        if struct.restraints is not None:
            if len(struct.restraints.rests) > 0:
                reslist = copy.deepcopy(struct.restraints)
                if restraints is not None:
                    for b in restraints:
                        found=False
                        for a in reslist.rests:
                            if a.is_same(b):
                                a.value=b.value
                                found=True
                        if not found:
                            reslist.rests.append(b)
        
        if reslist is None and restraints is not None:
            reslist = RestraintList( copy.deepcopy(restraints) )


        origatoms = struct.GetASEAtoms()
        myatoms = struct.GetASEAtoms()
        try:
            calc = get_persistent_calc(los, struct, reslist=reslist)
        except Exception:
            myatoms.calc = los.BuildRestrainedCalc(struct, reslist=reslist)
            calc = myatoms.calc
        else:
            myatoms.calc = calc

        conslist = None
        if struct.constraints is not None:
            if len(struct.constraints) > 0:
                conslist = copy.deepcopy(struct.constraints)
                if constraints is not None:
                    for b in constraints:
                        found=False
                        for a in conslist.cons:
                            if a.is_same(b):
                                a.value = b.value
                                found=True
                        if not found:
                            conslist.cons.append(b)
                            
        if conslist is None and constraints is not None:
            conslist = ConstraintList( copy.deepcopy(constraints) )
            
        asecons = None
        cons = None
        if conslist is not None:
            cons = conslist.FillConstraints(myatoms,force=False)
            origcons = conslist.FillConstraints(myatoms,force=True)
            myatoms = ApplyConstraints(myatoms,cons,graph=struct.GetGraph()) #,rests=reslist.rests,k=1)
            if cons is not None:
                asecons = to_ase(cons)


            
    
        del myatoms.constraints
        myatoms.set_constraint( asecons )
        myatoms.calc = calc
        myatoms.calc.reset()
        optimizer = BFGS(myatoms, logfile=None)
        optimizer.run(fmax=los.args.ase_opt_tol,steps=los.args.geometric_maxiter)
        
        ene = myatoms.get_potential_energy()
        crd = myatoms.get_positions()
        frc = myatoms.get_forces()
        out = copy.deepcopy(struct)
        out.Update(ene,crd,frc)
                                           

        if True:
            from . Constraints import FillConstraints
            if cons is not None:
                cvals = FillConstraints(out,cons,force=True)
                ovals = FillConstraints(origatoms,cons,force=True)
                for i in range(len(cvals)):
                    if cons[i].value is not None and \
                       cvals[i].value is not None and \
                       ovals[i].value is not None:
                        print("Constraint %2i tgt=%9.2f opt=%9.2f orig=%9.2f"%\
                              ( i+1, cons[i].value, cvals[i].value,
                                ovals[i].value ) )
            
        if True:
            if reslist is not None:
                ocrd = origatoms.get_positions()
                for i in range(len(reslist)):
                    val = reslist[i].GetCrdValue(crd)
                    oval = reslist[i].GetCrdValue(ocrd)
                    print("Restraint  %2i tgt=%9.2f obs=%9.2f orig=%9.2f"%\
                          ( i+1, reslist[i].value, val, oval ) )
            

        if cons is not None:
            out.constraints = ConstraintList( cons )
            out.data["constraints"] = out.constraints.to_list_of_dict()
        if reslist is not None:
            out.restraints = reslist
            out.data["restraints"] = out.restraints.to_list_of_dict()
    return out




def _run_geometric_with_watchdog(cmds, tmplog,
                                 poll_interval_sec=5.0,
                                 stall_timeout_sec=30 * 60,
                                 bmatrix_wedge_pattern="more than 1000 B-matrices stored"):
    """
    Execute geomeTRIC as an isolated background subprocess while monitoring log output logs.

    This routine protects optimization loops against internal coordinate parameter 
    bottlenecks or matrix calculation stalls. It polls log files for specific runaway 
    patterns or extended file lock freezes, issuing SIGTERM and SIGKILL sequences to 
    the active child process group if an anomaly occurs.

    Parameters
    ----------
    cmds : list of str
        The command sequence array mapping executable names and parameter option flags.
    tmplog : str
        The path tracking the active output stream file logs.
    poll_interval_sec : float, default 5.0
        The sleep period interval tracking how frequently file system sizes are pulled.
    stall_timeout_sec : float, default 1800.0
        The execution time limit tolerance checking for log file size stagnation.
    bmatrix_wedge_pattern : str, default 'more than 1000 B-matrices stored'
        The string warning pattern triggering immediate process group termination.

    Raises
    ------
    RuntimeError
        If matrix warnings are detected or if file system changes freeze over limits.
    """
    import os
    import signal
    import subprocess as subp
    import time

    child_env = os.environ.copy()
    extra = "ignore:.*ignore_bad_restart_file:FutureWarning"
    existing = child_env.get("PYTHONWARNINGS", "")
    if "ignore_bad_restart_file" not in existing:
        child_env["PYTHONWARNINGS"] = (
            f"{existing},{extra}" if existing else extra
        )


    
    proc = subp.Popen(cmds, text=True, env=child_env,
                      start_new_session=True)

    log_pos = 0
    # carry the tail of the last chunk so a pattern split across two
    # reads is still detected
    log_tail = ""
    last_change = time.monotonic()
    wedge_reason = None

    try:
        while proc.poll() is None:
            try:
                cur_size = os.path.getsize(tmplog)
            except OSError:
                cur_size = 0

            if cur_size > log_pos:
                try:
                    with open(tmplog, "r") as fh:
                        fh.seek(log_pos)
                        chunk = fh.read()
                except OSError:
                    chunk = ""
                if chunk:
                    log_pos += len(chunk)
                    scan_text = log_tail + chunk
                    if bmatrix_wedge_pattern in scan_text:
                        wedge_reason = (
                            f"geomeTRIC wedged: '{bmatrix_wedge_pattern}'"
                            f" detected in {tmplog}"
                        )
                        break
                    log_tail = scan_text[-len(bmatrix_wedge_pattern):]
                    last_change = time.monotonic()
            elif cur_size < log_pos:
                # log was rotated/truncated
                log_pos = 0
                log_tail = ""
                last_change = time.monotonic()
            elif time.monotonic() - last_change > stall_timeout_sec:
                wedge_reason = (
                    f"geomeTRIC stalled: no log growth in {tmplog} "
                    f"for {stall_timeout_sec}s"
                )
                break

            time.sleep(poll_interval_sec)
    finally:
        if proc.poll() is None:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            except (OSError, ProcessLookupError):
                pass
            try:
                proc.wait(timeout=10)
            except subp.TimeoutExpired:
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except (OSError, ProcessLookupError):
                    pass
                try:
                    proc.wait(timeout=5)
                except subp.TimeoutExpired:
                    pass

    if wedge_reason is not None:
        raise RuntimeError(wedge_reason)


def GeomOpt_GEOMETRIC(los,struct,constraints=None,restraints=None):
    """
    Perform a structural relaxation pass using the geomeTRIC optimization package.

    This function exports coordinates, writes the same constraint file as the
    CLI path, and by default runs geomeTRIC **in-process** against a persistent
    ASE calculator (avoids spawning a new interpreter and reloading the ML
    model at every node). Set ``FFPOPT_GEOMETRIC_SUBPROCESS=1`` to restore the
    legacy ``geometric-optimize`` subprocess + watchdog path.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The container tracking execution arguments and environmental parameters.
    struct : ffpopt.Struct.Struct
        The initial structure to track, sanitize, and minimize.
    constraints : list of Constraint, optional
        Active parameters tracking internal coordinate constraints to apply.
    restraints : list of Restraint, optional
        Active parameters tracking system restraint metrics to apply.

    Returns
    -------
    ffpopt.Struct.Struct
        The optimized molecular object with updated spatial coordinates.
    """

    import os
    from . CpuThreads import pin_math_threads
    pin_math_threads(1)

    import copy
    import ase.io
    from . constants import AU_PER_ELECTRON_VOLT
    from . Options import argparse2geometric, GetStandardOptions, configure_geometric_logging
    from . Constraints import ConstraintList
    from . Restraints import RestraintList
    from . Constraints import ApplyConstraints
    from . Struct import ListOfStruct
    from . AseEngine import (
        get_persistent_calc,
        run_geometric_inprocess,
        use_geometric_subprocess,
        prepare_geometric_tmpdir,
    )
    from tempfile import mkstemp
    import os
    import sys
    from pathlib import Path
    import subprocess as subp

    
    if True:

        tmpfile_loc = "./tmpfiles"
        if not Path(tmpfile_loc).is_dir():
            os.makedirs(tmpfile_loc, exist_ok=True)
        
        fd,tmpxyz = mkstemp(dir=tmpfile_loc,prefix="tmp.",suffix=".xyz")
        if not os.isatty(fd):  # Check if fd is still valid
            os.close(fd)

        tmpbase = str(Path(tmpxyz).with_suffix(""))
        tmpopt  = tmpbase + "_optim.xyz"
        tmplog  = tmpbase + ".log"
        tmpcons = tmpbase + ".cons.inp"
        tmpdir  = tmpbase + ".tmp"
        tmprst  = tmpbase + ".rst.inp"
        tmpjson = tmpbase + ".json"

        origatoms = struct.GetASEAtoms()
        myatoms = copy.deepcopy(origatoms)

        reslist = None
        if struct.restraints is not None:
            if len(struct.restraints.rests) > 0:
                reslist = copy.deepcopy(struct.restraints)
                if restraints is not None:
                    for b in restraints:
                        found=False
                        for a in reslist.rests:
                            if a.is_same(b):
                                a.value=b.value
                                found=True
                        if not found:
                            reslist.rests.append(b)
        
        if reslist is None and restraints is not None:
            reslist = RestraintList( copy.deepcopy(restraints) )
        
        conslist = None
        if struct.constraints is not None:
            if len(struct.constraints) > 0:
                conslist = copy.deepcopy(struct.constraints)
                
                if constraints is not None:
                    for b in constraints:
                        found=False
                        for a in conslist.cons:
                            if a.is_same(b):
                                a.value = b.value
                                found=True
                        if not found:
                            conslist.cons.append(b)
                            
        if conslist is None and constraints is not None:
            conslist = ConstraintList( copy.deepcopy(constraints) )


        cons = None
        if conslist is not None:
            cons = conslist.FillConstraints(myatoms,force=False)
            origcons = conslist.FillConstraints(myatoms,force=True)
            myatoms = ApplyConstraints(myatoms,cons,graph=struct.GetGraph()) #,rests=reslist.rests)


        if conslist is not None:
            fh = open(tmpcons,"w")
            clines = conslist.to_geometric()
            freezes = []
            nonfreezes = []
            for line in clines:
                if "xyz" in line[0:4]:
                    freezes.append(line)
                else:
                    nonfreezes.append(line)
            if len(nonfreezes) > 0:
                fh.write("$set\n")
                for line in nonfreezes:
                    fh.write("%s\n"%(line))
            if len(freezes) > 0:
                fh.write("$freeze\n")
                for line in freezes:
                    fh.write("%s\n"%(line))
            fh.close()

        inproc = not use_geometric_subprocess()
        if inproc:
            try:
                from geometric.ase_engine import EngineASE  # noqa: F401
            except ImportError:
                inproc = False

        crd = None
        ene = None
        frc = None
        if inproc:
            prepare_geometric_tmpdir(tmpbase)
            calc = get_persistent_calc(los, struct, reslist=reslist)
            myatoms.calc = calc
            try:
                calc.reset()
            except Exception:
                pass
            geo = GetStandardOptions(los.args).get("geometric", {}) or {}
            log_ini = None
            try:
                log_ini = configure_geometric_logging(
                    getattr(los.args, "geometric_ini", None)
                )
            except FileNotFoundError:
                log_ini = None
            if getattr(los.args, "geometric_ini", None) is not None:
                if len(str(los.args.geometric_ini)) == 0:
                    log_ini = ""
            result = run_geometric_inprocess(
                myatoms,
                calc,
                prefix=tmpbase,
                constraints_path=tmpcons if conslist is not None else None,
                coordsys=geo.get("coordsys", getattr(los.args, "geometric_coordsys", "tric")),
                maxiter=int(geo.get("maxiter", getattr(los.args, "geometric_maxiter", 500))),
                converge=geo.get("converge", getattr(los.args, "geometric_converge", "set GAU")),
                enforce=float(geo.get("enforce", getattr(los.args, "geometric_enforce", 0.0))),
                log_ini=log_ini if log_ini else None,
            )
            crd = result["coords"]
            if result["energy_ha"] is not None:
                ene = result["energy_ha"] / AU_PER_ELECTRON_VOLT()
            else:
                myatoms.set_positions(crd)
                myatoms.calc = calc
                ene = myatoms.get_potential_energy()
        else:
            mystruct = copy.deepcopy(struct)
            if conslist is not None:
                mystruct.constraints = None
                mystruct.data["constraints"] = []

            if reslist is not None:
                mystruct.restraints = reslist
                mystruct.data["restraints"] = reslist.to_list_of_dict()

            mylos = ListOfStruct( [mystruct] )
            mylos.save(tmpjson)

            cmds = argparse2geometric(tmpjson,los.args)
            cmds.append( tmpxyz )
            if conslist is not None:
                cmds.append(tmpcons)

            ase.io.write(tmpxyz,myatoms, format='xyz', parallel=False)

            _run_geometric_with_watchdog(cmds, tmplog)
            if not os.path.exists(tmpopt):
                raise Exception(f"File not found: {tmpopt}")

            out_atoms = ase.io.read(tmpopt,index='-1',parallel=False)
            out_atoms.set_initial_charges( myatoms.get_initial_charges() )
            keys = [ key for key in out_atoms.info ]
            ene = float(keys[-1]) / AU_PER_ELECTRON_VOLT()
            crd = out_atoms.get_positions()

        out = copy.deepcopy(struct)
        out.Update(ene,crd,frc)


        if True:
            from . Constraints import FillConstraints
            if cons is not None:
                cvals = FillConstraints(out,cons,force=True)
                ovals = FillConstraints(origatoms,cons,force=True)
                for i in range(len(cvals)):
                    if cons[i] is not None and \
                       cvals[i].value is not None and \
                       ovals[i] is not None:
                        print("Constraint %2i tgt=%9.2f opt=%9.2f orig=%9.2f"%\
                              ( i+1, cons[i].value, cvals[i].value,
                                ovals[i].value ) )
            
        if True:
            if reslist is not None:
                ocrd = origatoms.get_positions()
                for i in range(len(reslist)):
                    val = reslist[i].GetCrdValue(crd)
                    oval = reslist[i].GetCrdValue(ocrd)
                    print("Restraint  %2i tgt=%9.2f obs=%9.2f orig=%9.2f"%\
                          ( i+1, reslist[i].value, val, oval ) )
            
        
        if cons is not None:
            out.constraints = ConstraintList( cons )
            out.data["constraints"] = out.constraints.to_list_of_dict()

        if reslist is not None:
            out.restraints = reslist
            out.data["restraints"] = out.restraints.to_list_of_dict()
            
        for f in [tmpxyz,tmpopt,tmplog,tmpcons]:
            if os.path.exists(f):
                os.remove(f)

        if os.path.isdir(tmpdir):
            import shutil
            shutil.rmtree(tmpdir)
    
    return out



def GeomOpt_SinglePoint(los,struct,constraints=None,restraints=None):
    """
    Perform a single-point electronic energy evaluation without relaxing positions.

    This function sets up a single potential evaluation pass on the unrelaxed input coordinates,
    updates constraint and restraint definitions, and returns the structural instance loaded
    with accurate energy and gradient force tensors.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The orchestrating batch manager object tracking operational variables.
    struct : ffpopt.Struct.Struct
        The target molecular geometry object to query.
    constraints : list of Constraint, optional
        A collection of configuration parameters tracking coordinate boundaries.
    restraints : list of Restraint, optional
        A collection of parameter configurations tracking restraints to evaluate.

    Returns
    -------
    ffpopt.Struct.Struct
        The output structure equipped with updated energy and force parameters.
    """
    import os
    import copy
    from . Constraints import ConstraintList
    from . Constraints import to_ase,ApplyConstraints
    from . Restraints import RestraintList

    if True:
        reslist = None
        if struct.restraints is not None:
            if len(struct.restraints.rests) > 0:
                reslist = copy.deepcopy(struct.restraints)
                if restraints is not None:
                    for b in restraints:
                        found=False
                        for a in reslist.rests:
                            if a.is_same(b):
                                a.value=b.value
                                found=True
                        if not found:
                            reslist.rests.append(b)
        
        if reslist is None and restraints is not None:
            reslist = RestraintList( copy.deepcopy(restraints) )
        
        
        myatoms = struct.GetASEAtoms()
        try:
            from . AseEngine import get_persistent_calc
            myatoms.calc = get_persistent_calc(los, struct, reslist=reslist)
        except Exception:
            myatoms.calc = los.BuildRestrainedCalc(struct,reslist=reslist)
        calc = myatoms.calc

        conslist = None
        if struct.constraints is not None:
            if len(struct.constraints) > 0:
                conslist = copy.deepcopy(struct.constraints)
                if constraints is not None:
                    for b in constraints:
                        found=False
                        for a in conslist.cons:
                            if a.is_same(b):
                                a.value = b.value
                                found=True
                        if not found:
                            conslist.cons.append(b)
        if conslist is None and constraints is not None:
            conslist = ConstraintList( copy.deepcopy(constraints) )

        asecons = None
        cons = None
        if conslist is not None:
            cons = conslist.FillConstraints(myatoms)
            myatoms = ApplyConstraints(myatoms,cons,graph=struct.GetGraph())
            asecons = to_ase(cons)

        del myatoms.constraints
        myatoms.set_constraint( asecons )
        myatoms.calc = calc
        myatoms.calc.reset()

        ene = myatoms.get_potential_energy()
        crd = myatoms.get_positions()
        frc = myatoms.get_forces()
        out = copy.deepcopy(struct)
        out.Update(ene,crd,frc)
        
        if cons is not None:
            out.constraints = ConstraintList( cons )
            out.data["constraints"] = out.constraints.to_list_of_dict()
        if reslist is not None:
            out.restraints = reslist
            out.data["restraints"] = out.restraints.to_list_of_dict()
            
    return out



def GeomOpt(los,struct,constraints=None,restraints=None):
    """
    Direct and execute a calculation task matching the selected CLI option flags.

    This routing function acts as a centralized dispatcher. It reads the command line 
    arguments embedded inside the list of structures to choose between a basic 
    single-point calculation, a standard geomeTRIC optimization run, or a primary 
    ASE BFGS relaxation pass equipped with an automatic geomeTRIC exception fallback path.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The container tracking runtime arguments, options, and environmental variables.
    struct : ffpopt.Struct.Struct
        The initial input molecular geometry and parameter context to evaluate.
    constraints : list of Constraint, optional
        A sequence of geometric constraints to merge and apply during execution.
    restraints : list of Restraint, optional
        A sequence of structural restraints to merge and apply during execution.

    Returns
    -------
    ffpopt.Struct.Struct
        The resulting molecular structure loaded with updated energy and force arrays.
    """

    if los.args.no_opt:
        out = GeomOpt_SinglePoint(los,struct,constraints,restraints)
    elif not los.args.geometric_opt:
        try:
            out = GeomOpt_ASE(los,struct,constraints,restraints)
        except Exception as e:
            import traceback
            print("\n\n\nASE GEOMETRY OPTIMIZATION FAILURE\n")
            print(e)
            traceback.print_exc()
            out = GeomOpt_GEOMETRIC(los,struct,constraints,restraints)
    else:
        out = GeomOpt_GEOMETRIC(los,struct,constraints,restraints)
    return out





def CheckForces(los,struct,delta=1.e-2):
    """
    Compare analytical force tensors against numerical finite-difference gradients.

    This verification driver performs a central finite-difference displacement 
    pass across all spatial dimensions for every atomic center. It evaluates 
    single-point energies at perturbed geometries, computes numerical forces, 
    and prints a detailed side-by-side comparison report to standard output.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The container tracking runtime options, calculators, and system definitions.
    struct : ffpopt.Struct.Struct
        The molecular configuration state instance undergoing gradient validation.
    delta : float, default 1e-2
        The displacement coordinate step step length expressed in Angstroms.

    Returns
    -------
    None
    """
    import copy
    import numpy as np
    from . constants import AU_PER_KCAL_PER_MOL, AU_PER_ELECTRON_VOLT

    DEL = delta
    ana = GeomOpt_SinglePoint(los,struct)
    fana = ana.get_forces()
    fnum = np.zeros( fana.shape )
    ts = copy.deepcopy(struct)
    n = len(struct.data["elements"])
    
    for a in range(n):
        for k in range(3):
            ts.data["positions"][a][k] += DEL
            o = GeomOpt_SinglePoint(los,ts)
            ehi = o.get_potential_energy()
            ts.data["positions"][a][k] -= 2*DEL
            o = GeomOpt_SinglePoint(los,ts)
            elo = o.get_potential_energy()
            ts.data["positions"][a][k] += DEL
            fnum[a,k] = - (ehi-elo)/(2*DEL)

    print("%4s %35s  %35s  %35s"%("Atom","Analytic Frc","Numerical Frc","Ana-Num"))
    for a in range(n):
        print("%5i %11.2e %11.2e %11.2e  %11.2e %11.2e %11.2e  %11.2e %11.2e %11.2e"%\
              ( a+1,
                fana[a,0],fana[a,1],fana[a,2],
                fnum[a,0],fnum[a,1],fnum[a,2],
                fana[a,0]-fnum[a,0],fana[a,1]-fnum[a,1],fana[a,2]-fnum[a,2] ))

            




def ApplyDihedConstraint(atoms,idxs,value,rotmask):
    """ Apply a dihedral constraint to an ASE Atoms object.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The geometry to modify.
    idxs : list of int
        A list of four 0-based atom indices defining the dihedral.
    value : float
        The dihedral angle in degrees to set.
    rotmask : list of bool
        A list of booleans indicating which atoms to rotate.
    
    Returns
    -------
    ase.Atoms
        A new ASE Atoms object with the modified dihedral angle.
    
    """
    import copy
    out = atoms.copy()
    out.set_dihedral(idxs[0],idxs[1],idxs[2],idxs[3],value,
                     mask=rotmask)
    return out


    

def DihedScan(los,struct,con,sched):
    """ Perform a dihedral scan using forward and reverse scans to find the minimum energy conformation.
    
    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        Structure list object used to build calculator
    struct : ffpopt.Struct.Struct
        A Struct object (usually los.structs[0])
    con : Constraint
        The dihedral constraint to scan.
    sched : list of float
        A list of dihedral angles in degrees to scan.
    
    Returns
    -------
    ListOfStruct
        A list of Structs corresponding to the scanned geometries. The Struct names
        are angXXX, where XXX is a zero-padded integer of the dihedral angle in degrees.
    """
    import numpy as np
    import copy
    from . Struct import ListOfStruct
    
    idxs = copy.deepcopy( con.idxs )

    mingeom = GeomOpt(los,struct)

    sched = np.array(sched,copy=True)
    minang = mingeom.get_dihedral(idxs[0],idxs[1],idxs[2],idxs[3])

    difangs = [ abs((x-minang)%360) for x in sched ]
    istart = np.argmin(difangs)
    nscan = len(sched)
    scan = []

    for i in range(nscan+1):
        icur = (istart + i) % nscan
        if icur == nscan:
            value = sched[0]
        else:
            value = sched[icur]

        if i == 0:
            igeom = mingeom.copy()
        else:
            igeom = opt.copy()

        fang = igeom.get_dihedral(idxs[0],idxs[1],idxs[2],idxs[3])
            
        cons = [ copy.deepcopy(con) ]
        cons[0].value = value

        opt = GeomOpt(los,igeom,constraints=cons)
        

        print("Finished angle:",value)

        opt.data["name"] = "ang%03i"%(int(round(value)))
        if i == nscan:
            if opt.get_potential_energy() < scan[0].get_potential_energy():
                scan[0] = opt
        else:
            scan.append(opt)

    scan = sorted(scan,key=lambda x: x.data["name"])

    return ListOfStruct(scan)


def FwdRevDihedScan_worker(args):
    """ Worker function for parallel forward and reverse dihedral scans."""
    return DihedScan(*args)


def FwdRevDihedScan(los,struct,con,sched,parallel=False):
    """ Perform a dihedral scan using both forward and reverse scans to find the minimum energy conformation.
    
    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        Structure list object used to build calculator
    struct : ffpopt.Struct.Struct
        A Struct object (usually los.structs[0])
    con : ffpopt.Constraint.Constraint
        The dihedral constraint to scan.
    sched : list of float
        A list of dihedral angles in degrees to scan.
    parallel : bool, optional
        Whether to perform the forward and reverse scans in parallel. Default is False.
    
    Returns
    -------
    ListOfStruct
        A list of Structs corresponding to the scanned geometries. The Struct names
        are angXXX, where XXX is a zero-padded integer of the dihedral angle in degrees.
    """
    
    import ase.io
    from . Struct import ListOfStruct
    
    revsched = sched[::-1]

    if parallel:
        proclist = [ (los,struct,con,sched),
                     (los,struct,con,revsched) ]

        from . NondaemonPool import make_nondaemon_spawn_pool
        pool = make_nondaemon_spawn_pool(2)
        try:
            olists = pool.map(FwdRevDihedScan_worker,proclist)
        finally:
            pool.close()
            pool.join()
        fwd = olists[0]
        rev = olists[1]
    else:
        fwd = DihedScan(los,struct,con,sched)
        rev = DihedScan(los,struct,con,revsched)

        
    scan = []
    
    for a,b in zip(fwd,rev):
        if a.data["name"] != b.data["name"]:
            raise Exception("Expected fwd and rev scans at the same set of dihedrals, but found %s and %s\n"%(a.data["name"],b.data["name"]))
        if a.data["energy"] < b.data["energy"]:
            scan.append(a)
        else:
            scan.append(b)

    return ListOfStruct(scan)




###########################################################################################
###########################################################################################
###########################################################################################
        

def ParallelGeomOpt(los,norestene,nproc):
    """
    Orchestrate molecular geometry optimizations across parallel worker environments.

    This router inspects the active execution architecture, automatically dispatching 
    the list of structures to an MPI-based cluster loop if an active communicator 
    context is detected, or falling back to a local multi-threaded ProcessPoolExecutor 
    sandbox topology layout.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The container bundle tracking the input configurations and evaluation arguments.
    norestene : bool
        If True, strips constraints and restraints from the final evaluation pass.
    nproc : int
        The maximum number of concurrent processor worker threads to spawn locally.

    Returns
    -------
    ffpopt.Struct.ListOfStruct
        A new collection holding the fully relaxed and evaluated molecular frames.
    """
    out = None
    if is_mpi():
        out = ParallelGeomOpt_mpi(los, norestene)
    else:
        out = ParallelGeomOpt_threads(los,norestene,nproc)
    return out
        

###########################################################################################
###########################################################################################
###########################################################################################
        
class CalcNode(object):
    """
    Isolated execution tracking container mapping to a single evaluation task.

    This container packages a single molecular structure frame alongside its 
    associated execution flags and reference calculators, serving as a unified 
    picklable vehicle passed safely down multiprocessing or MPI communication tracks.

    Attributes
    ----------
    los : ffpopt.Struct.ListOfStruct
        The base options and argument context used to build structural calculators.
    s : ffpopt.Struct.Struct
        The initial input geometry tracking state assigned to this worker node.
    norestene : bool
        If True, strips out active constraints and restraints during final steps.
    out : ffpopt.Struct.Struct or None
        The fully updated optimization result instance state, or None if uncalculated.
    """
    def __init__(self,los,s,norestene):
        """
        Initialize a calculation worker node tracking template instance.

        Parameters
        ----------
        los : ffpopt.Struct.ListOfStruct
            The argument context used to assign structural calculators.
        s : ffpopt.Struct.Struct
            The target molecular configuration state bound to this node.
        norestene : bool
            Operational parameter determining constraint stripping modes.
        """
        self.los = los
        self.s = s
        self.norestene = norestene
        self.out = None

    def calculate(self):
        """
        Execute the structural optimization and single-point property cleanup.

        Runs the core minimization pass using the designated engine flags. If 
        `norestene` is enabled, constraints and restraints are stripped from a 
        duplicate instance frame right before triggering a final single-point 
        evaluation to calculate clean, unconstrained target potential energies 
        and force matrix parameters.

        Returns
        -------
        None
        """
        from ffpopt.GeomOpt import GeomOpt,GeomOpt_SinglePoint
        import copy
        if self.los.args.no_opt:
            self.out = copy.deepcopy(self.s)
        else:
            self.out = GeomOpt(self.los,self.s)
        tmp = copy.deepcopy(self.out)
        if self.norestene:
            tmp.restraints = None
            tmp.constraints = None
        tmp = GeomOpt_SinglePoint(self.los,tmp)
        self.out.Update( tmp.get_potential_energy(), tmp.get_positions(), tmp.get_forces() )

        
def _run_node( node ):
    """
    Top-level pickle-compliant wrapper function to execute a worker node task.

    Parameters
    ----------
    node : CalcNode
        The input structural calculation node instance passed from executors.

    Returns
    -------
    CalcNode
        The updated worker node packed with completed results.
    """
    node.calculate()
    return node


def is_mpi_worker():
    """
    Query the MPI environment to check if the current process is an active rank worker.

    Returns
    -------
    bool
        True if the active world size is greater than 1 and rank is greater than 0.
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    return size > 1 and rank > 0


def is_mpi():
    """
    Query the environment to check if an active MPI cluster context is deployed.

    Returns
    -------
    bool
        True if the global world communicator size tracking variable is greater than 1.
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    #rank = comm.Get_rank()
    size = comm.Get_size()
    return size > 1


def ParallelGeomOpt_threads(los,norestene,nproc):
    """
    Execute optimization passes across local worker threads using ProcessPoolExecutor.

    This function instantiates independent `CalcNode` object containers for every 
    packed structure, distribution-maps them across local worker pools, and 
    reassembles the finished nodes into a new structural collection matrix list.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The compilation collection tracking molecular geometry target frames.
    norestene : bool
        Operational flag tracking whether parameters strip restraints at exit points.
    nproc : int
        The maximum number of local worker processor channels to instantiate.

    Returns
    -------
    ffpopt.Struct.ListOfStruct
        A newly instantiated collection holding the completed, optimized structures.
    """
    import concurrent.futures
    import multiprocessing
    from . Struct import ListOfStruct
    from . NondaemonPool import make_nondaemon_spawn_pool

    nodes = [ CalcNode(los,s,norestene) for s in los ]
    pool = make_nondaemon_spawn_pool(max(1, int(nproc)))
    try:
        results = pool.map(_run_node, nodes)
    finally:
        pool.close()
        pool.join()
    return ListOfStruct( [ node.out for node in results ] )


# -------------------------------------------------------------------------
# WORKER SIDE CAR ENVIRONMENT
# -------------------------------------------------------------------------
# Worker-level global storage variables
_WORKER_LOS = None
_WORKER_NORESTENE = None

def _worker_init(los, norestene):
    """
    Initialize background global state trackers on a joining cluster node.

    This routine caches configuration parameters and base list of structures 
    instances into the worker process's isolated global data layer. This step 
    enables subsequent structural iterations to map against pre-allocated variables 
    and prevents redundant network data synchronization steps across processing ranks.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The base configuration parameters template cached into worker memory space.
    norestene : bool
        The active operational flag determining constraint stripping behavior.

    Returns
    -------
    None
    """
    global _WORKER_LOS, _WORKER_NORESTENE
    _WORKER_LOS = los
    _WORKER_NORESTENE = norestene

def _run_node_mpi(s):
    """
    Execute a single structural minimization using cached background metadata.

    This routine unpacks a standalone structure state configuration instance passed 
    over the network mesh network, instantiates a local execution CalcNode tracker, 
    triggers the optimization routines, and extracts the finalized structural 
    payload to minimize MPI data footprint overhead.

    Parameters
    ----------
    s : ffpopt.Struct.Struct
        The raw coordinate layout structure state instance sent by the master node.

    Returns
    -------
    ffpopt.Struct.Struct
        The localized execution output structure loaded with updated energies and forces.
    """
    global _WORKER_LOS, _WORKER_NORESTENE
    
    # Instantiate the node locally using the cached background variables
    node = CalcNode(_WORKER_LOS, s, _WORKER_NORESTENE)
    node.calculate()
    
    # Return ONLY the structure output payload to minimize MPI data footprint
    return node.out

# -------------------------------------------------------------------------
# TARGET MPI FUNCTION
# -------------------------------------------------------------------------
def ParallelGeomOpt_mpi(los, norestene):
    """
    Execute batch optimizations across a distributed cluster using an MPI worker pool.

    This driver uses MPICommExecutor to partition the COMM_WORLD space. Rank 0 assumes 
    the role of master task scheduler, manually broadcasting initialization payloads to 
    set up global worker variables before streaming individual structural configurations 
    across the active cluster mesh network with perfect load-balancing.

    Parameters
    ----------
    los : ffpopt.Struct.ListOfStruct
        The batch structural collection dataset containing targets to distribute.
    norestene : bool
        Operational flag tracking whether exit coordinates strip constraints.

    Returns
    -------
    ffpopt.Struct.ListOfStruct or None
        A new structural collection holding relaxed results on Rank 0, or None on 
        passive worker ranks.
    """
    from mpi4py import MPI
    from mpi4py.futures import MPICommExecutor
    from . Struct import ListOfStruct
    
    # MPICommExecutor partitions COMM_WORLD.
    # Workers enter a passive processing loop inside the 'with' block context.
    # Only Rank 0 exits the block to submit jobs via the executor.
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            # Set up global contextual environments on worker memory pools
            # Note: MPICommExecutor does not support the 'initializer' parameter, 
            # so we map the initialization function across workers manually.
            num_workers = MPI.COMM_WORLD.Get_size() - 1
            if num_workers > 0:
                list(executor.map(_worker_init, [los]*num_workers, [norestene]*num_workers))

            # Dynamically stream data chunks to achieve perfect load balancing
            results_iterator = executor.map(_run_node_mpi, list(los), chunksize=1)
            final_outputs = list(results_iterator)
            out = ListOfStruct(final_outputs)
            return out
    return None
