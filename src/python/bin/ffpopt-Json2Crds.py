#!/usr/bin/env python3


if __name__ == "__main__":

    
    from pathlib import Path
    import numpy as np
    import argparse
    import json
    
    import parmed as pmd
    from parmed.structure import Structure
    from parmed.topologyobjects import Atom, Bond, Residue
    
    import ase
    import ase.io
    from ffpopt.Struct import ListOfStruct
    from ffpopt.constants import GetAtomicNumber
    
    parser = argparse.ArgumentParser \
        ( formatter_class=argparse.RawDescriptionHelpFormatter,
          description="""Write coordinates""" )

    
    parser.add_argument \
        ("--inp","-i",
         help="input json file",
         required=True,
         type=str)

    parser.add_argument \
        ("--name",
         help="Only extract structure with this name",
         type=str)
    
    parser.add_argument \
        ("--out","-o",
         help="""output crd file. Suffix '.xyz' '.mol2' '.rst7' corresponds to xyz, mol2, and amber restart.
         tdep-inchi:        tautomer-dependent inchi string (nonstandard inchi)
         tindep-inchi:      tautomer-independent inchi string (standard inchi)
         tdep-smiles:       tautomer-dependent smiles string
         tindep-smiles:     tautomer-independent smiles string
         tdep-inchikey:     tautomer-dependent inchikey (only output)
         tindep-inchikey:   tautomer-independent inchikey (only output)
         tdep-inchitoken:   url-safe base64 tdep-inchi
         tindep-inchitoken: url-safe base64 tindep-inchi
         iupac: iupac molecule name""",
         required=True,
         type=str)
    
    args = parser.parse_args()

    opath = Path(args.out)
    
    inps = ListOfStruct.from_file(args.inp)
    found=False
    names = []
    for inp in inps:
        names.append(inp.data["name"])
        if args.name is not None:
            if inp.data["name"] != args.name:
                continue
            else:
                found=True
        if len(inps.structs) > 1 and args.name is None:
            base = opath.with_suffix('')
            s = opath.suffix
            s = "_" + inp.data["name"] + str(s)
            o = str(base) + s
        else:
            o = str(opath)
        if args.out.lower() == "tdep-inchi":
            print( inp.GetInchi(tautdep=True) )
        elif args.out.lower() == "tindep-inchi":
            print( inp.GetInchi(tautdep=False) )
        elif args.out.lower() == "tdep-inchikey":
            print( inp.GetInchiKey(tautdep=True) )
        elif args.out.lower() == "tindep-inchikey":
            print( inp.GetInchiKey(tautdep=False) )
        elif args.out.lower() == "tdep-inchitoken":
            print( inp.GetInchiToken(tautdep=True) )
        elif args.out.lower() == "tindep-inchitoken":
            print( inp.GetInchiToken(tautdep=False) )
        elif args.out.lower() == "tdep-smiles":
            print( inp.GetSmiles(tautdep=True) )
        elif args.out.lower() == "tindep-smiles":
            print( inp.GetSmiles(tautdep=False) )
        elif args.out.lower() == "iupac":
            print( inp.GetIUPAC() )

        else:
            inp.SaveCrds(o)

    if args.name is not None and not found:
        raise Exception(f"Could not find structure with name {args.name} in {names}")
        
        # eles  = inp.data["elements"]
        # types = inp.data["types"]
        # qs    = inp.data["charges"]
        # q     = inp.data["charge"]
        # crds  = inp.data["positions"]
        # bonds = inp.data["bonds"]
        
        
        # if opath.suffix == ".mol2":
        #     mol = Structure()
        #     alist = []
        #     for i in range(len(eles)):
        #         z = GetAtomicNumber(eles[i])
        #         a = Atom(name=eles[i],type=types[i],
        #                  atomic_number=z,charge=qs[i])
        #         a.xx = crds[i][0]
        #         a.xy = crds[i][1]
        #         a.xz = crds[i][2]
        #         alist.append(a)
        #         mol.add_atom(a,resname='MOL', resnum=1)
        #     for x in bonds:
        #         mol.bonds.append(pmd.Bond(alist[x[0]],alist[x[1]]))
        #     base = str(opath).replace(".mol2","_%s"%(inp.data["name"]))
        #     fname = base + ".mol2"
        #     mol.save(fname,overwrite=True)
        # elif opath.suffix == ".xyz":
        #     base = str(opath).replace(".xyz","_%s"%(inp.data["name"]))
        #     fname = base + ".xyz"
        #     atoms = inp.GetASEAtoms()
        #     ase.io.write(fname,atoms,format="extxyz")
        # elif opath.suffix == ".rst7":
        #     base = str(opath).replace(".rst7","_%s"%(inp.data["name"]))
            
    #import os
    #os._exit(0)

    
