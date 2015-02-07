#!/usr/bin/env python
# encoding: utf-8

import Bio.PDB as bp
import sys
import numpy as np

PDB_file = sys.argv[1]

def hbond_angle_term( Nx, Ns_COx, Ns_CHx, Cx, Cs_CHx, Cs_NHx ):
    vec_CO_N = Nx - Ns_COx
    vec_CH_N = Nx - Ns_CHx

    vec_N_H = vec_CO_N/np.linalg.norm(vec_CO_N) + vec_CH_N/np.linalg.norm(vec_CH_N)
    vec_N_H = vec_N_H/np.linalg.norm(vec_N_H)

    vec_CH_C = Cx - Cs_CHx
    vec_NH_C = Cx - Cs_NHx

    vec_C_O = vec_CH_C/np.linalg.norm(vec_CH_C) + vec_NH_C/np.linalg.norm(vec_NH_C)
    vec_C_O = vec_C_O/np.linalg.norm(vec_C_O)

    vec_N_C = Cx - Nx
    vec_N_C = vec_N_C/np.linalg.norm(vec_N_C)

    costhetaN = np.dot( vec_N_H, vec_N_C )
    costhetaC = np.dot( vec_C_O, -vec_N_C )

    if costhetaC > 0.0 and costhetaN > 0.0:
        return costhetaC*costhetaC*costhetaN*costhetaN

    return 0.0

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

    #Build list of pairs qualifying on:
    #not pro, within dist cuttoff, not adjacent or same residue
    for N in Ns:

        #check N doesn't belong to Pro; can't Hbond
        if N.get_parent().get_resname() == '  P':
            continue

        for C in Cs:

            #Within cutoff?
            if np.linalg.norm( N.get_coord() - C.get_coord() ) > hb_dist_cut:
                continue

            #Not on adjacent or the same residues?
            dist = abs( N.parent.get_id()[1] - C.parent.get_id()[1] )
            same_chain = ( N.parent.parent.get_id() == C.parent.parent.get_id() )

            if same_chain and dist < 2:
                continue

            pairs_within_cutoff.append( [N,C] )

    for N,C in pairs_within_cutoff:
        Ns_res = N.parent
        Cs_res = C.parent

        Ns_chain = Ns_res.parent
        Cs_chain = Cs_res.parent

        N_resid=Ns_res.get_id()[1]
        C_resid=Cs_res.get_id()[1]

        try:
            Ns_CO = filter( lambda x: x.name == 'C', Ns_chain[N_resid-1] )[0]
            Cs_NH = filter( lambda x: x.name == 'N', Cs_chain[C_resid+1] )[0]
        except KeyError:
            continue

        Ns_CH = filter(lambda x : x.name == 'CA', Ns_res)[0]
        Cs_CH = filter(lambda x : x.name == 'CA', Cs_res)[0]

        angle_term = hbond_angle_term( N.get_coord(), Ns_CO.get_coord(), Ns_CH.get_coord(),
                          C.get_coord(), Cs_CH.get_coord(), Cs_NH.get_coord()  )
        dist = np.linalg.norm( N.get_coord() - C.get_coord() )
        dist_term = 5*(sigma/dist)**12 - 6*(sigma/dist)**10

        pre_energy = -dist_term*angle_term

        if pre_energy > 0.5:
            print "{0:2d} {1:2d}: {2:.2f}%".format(N_resid, C_resid, 100*pre_energy)

    return pairs_within_cutoff

if __name__ == '__main__':
    for i_model, model in enumerate(bp.PDBParser().get_structure("",PDB_file)):

        print "-------------------"
        print "Model", i_model
        print "-------------------"
        pairs =  find_hbonds( model )
