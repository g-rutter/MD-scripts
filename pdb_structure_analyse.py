#!/usr/bin/env python
# encoding: utf-8

import Bio.PDB as bp
import sys
import numpy as np

PDB_file = sys.argv[1]

def find_hbonds( model ):

    #settings
    hb_dist_cut=8.0
    sigma=4.11
    e_hb=3.60

    #Create list of NH-
    Ns = []
    Cs = []
    pairs_within_cutoff = []
    for atom in model.get_atoms():
        if atom.name == 'N':
            Ns.append(atom)
        elif atom.name == 'C':
            Cs.append(atom)

    #Build list of pairs qualifying on all criteria except angles.
    for i_N, N in enumerate(Ns):

        #check N doesn't belong to Pro; can't Hbond
        if N.get_parent().get_resname() == '  P':
            continue

        for C in Cs[i_N:]:

            #Within cutoff?
            if np.linalg.norm( N.get_coord() - C.get_coord() ) > hb_dist_cut:
                continue

            #Not on adjacent or the same residues?
            dist = abs( N.parent.get_id()[1] - C.parent.get_id()[1] )
            same_parent = ( N.parent.parent.get_id() == C.parent.parent.get_id() )

            if same_parent and dist < 2:
                continue

            pairs_within_cutoff.append( [N,C] )

    for N,C in pairs_within_cutoff:
        N_resid=N.parent.get_id()[1]
        C_resid=C.parent.get_id()

        try:
            print C_resid
            Cs_lN = C.parent.parent.get_list()[C_resid[1]+1]
            print Cs_lN
        except:
            pass



    return pairs_within_cutoff

if __name__ == '__main__':
    for i_model, model in enumerate(bp.PDBParser().get_structure("",PDB_file)):

        pairs =  find_hbonds( model )
