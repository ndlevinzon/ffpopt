#!/usr/bin/env python3

from collections import defaultdict as ddict

def parmed2ase(mol):
    from parmed import periodic_table
    import numpy as np
    import ase
    import numpy as np
    qs = np.array([ a.charge for a in mol.atoms ])
    qsum = sum(qs)
    charge = int(round(sum([ a.charge for a in mol.atoms ])))
    qs += (charge-qsum)/len(qs)
    eles = [ periodic_table.Element[a.element]
            for a in mol.atoms]
    crds = np.array([ [ a.xx, a.xy, a.xz ] for a in mol.atoms ])
    atlist = "".join( ["%s1"%(ele) for ele in eles ] )
    
    return ase.Atoms(atlist,positions=crds,charges=qs),charge
    


def parmed2graph(mol):
    import parmed
    import numpy as np
    from . GraphSearch import GraphSearch
    edges = []
    for x in mol.bonds:
        edges.append( "%i~%i"%(x.atom1.idx,x.atom2.idx) )
    return GraphSearch(edges)



def bonds2graph(bonds):
    from . GraphSearch import GraphSearch
    edges = []
    for x in bonds:
        edges.append( "%i~%i"%(x[0],x[1]) )
    return GraphSearch(edges)



def CopyParm( parm ):
    import copy
    try:
        parm.remake_parm()
    except:
        pass
    p = copy.copy( parm )
    p.coordinates = copy.copy( parm.coordinates )
    p.box = copy.copy( parm.box )
    try:
        p.hasbox = copy.copy( parm.hasbox )
    except:
        p.hasbox = False
    return p


def RotateMask(graph,idxs):
    left = "%i"%(idxs[1])
    right = "%i"%(idxs[2])
    nat = len(graph.nodes)
    gleft = []
    gright = []
    for n in graph.nodes:
        if n == left:
            gleft.append(int(n))
        elif n == right:
            gright.append(int(n))
        else:
            rleft = len(graph.FindMinPaths(left,n)[0])
            rright = len(graph.FindMinPaths(right,n)[0])
            if rleft < rright:
                gleft.append(int(n))
            else:
                gright.append(int(n))
    mask = [0]*nat
    if len(gleft) < len(gright):
        gmove = gleft
    else:
        gmove = gright
    for i in gmove:
        mask[i] = 1

    if mask[idxs[3]] == 0:
        mask = [ 1-x for x in mask ]
        
    return mask



def RotateBondMask(graph,bondpair):
    left = "%i"%(bondpair[0])
    right = "%i"%(bondpair[1])
    nat = len(graph.nodes)
    gleft = []
    gright = []
    for n in graph.nodes:
        if n == left:
            gleft.append(int(n))
        elif n == right:
            gright.append(int(n))
        else:
            rleft = len(graph.FindMinPaths(left,n)[0])
            rright = len(graph.FindMinPaths(right,n)[0])
            if rleft < rright:
                gleft.append(int(n))
            else:
                gright.append(int(n))
    mask = [0]*nat
    for i in gleft:
        mask[i] = 1
    
    return mask




def CheckParmCMAP(file_path):
    """
    Checks an Amber parm7 file for CMAP entries. If either value in 
    CMAP_COUNT is 0, duplicates the original file with a '.bak' extension, 
    removes the target CMAP sections, and overwrites the original file.
    """
    import parmed as pmd
    import shutil

    # Load as a raw Amber topology object to retain header metadata intact
    try:
        p = pmd.amber.AmberParm(file_path)
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return

    # Check if the target flag exists in the topology backend
    if "CMAP_COUNT" in p.parm_data:
        cmap_counts = p.parm_data["CMAP_COUNT"]
        
        # Check if either count inside the data block evaluates to 0
        if any(count == 0 for count in cmap_counts):
            
            # 1. Create a safe duplicate backup of the file before altering it
            backup_path = f"{file_path}.bak"
            try:
                shutil.copy2(file_path, backup_path)
                print(f"Created a secure backup copy at: {backup_path}")
            except Exception as backup_error:
                print(f"Aborting process: Failed to create backup file. Error: {backup_error}")
                return

            # 2. Define target metadata blocks to remove
            flags_to_remove = ["CMAP_COUNT", "CMAP_RESOLUTION", "CMAP_INDEX"]
            
            for flag in flags_to_remove:
                # Remove raw numerical arrays
                if flag in p.parm_data:
                    del p.parm_data[flag]
                # Remove corresponding print formats (%FORMAT)
                if flag in p.formats:
                    del p.formats[flag]
                # Clean up cached file comments if any exist
                #if flag in p.comments:
                #    del p.comments[flag]
                # Erase entry order index so the blocks are skipped entirely
                if flag in p.flag_list:
                    p.flag_list.remove(flag)
                    
            # 3. Overwrite the original structure safely
            p.write_parm(file_path)
            print(f"Altered topology written: Removed CMAP records from {file_path}.")
        #else:
        #    print("CMAP_COUNT fields are non-zero. No actions taken.")
    #else:
    #    print("No CMAP_COUNT section found inside this parm7 topology.")
