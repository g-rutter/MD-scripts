#!/usr/bin/env python
# encoding: utf-8

import Bio.PDB as bp
import sys
from numpy import pi, array, linalg, dot, zeros, set_printoptions
from collections import defaultdict

set_printoptions(linewidth=200)

def count_in_area(area, coordinates):
    """
    area is an area specified as ( ( x_low, x_high ), ( y_low, y_high ) )
    coordinates is an iterable containing N objects, each of which is an
    iterable containing two items forming a 2D coordinate.

    Returns the number of those coordinates which fall within the area.
    """

    npcoords = array(coordinates)

    #Check if each value is in its range
    #i.e. xlo < x < xhi and ylo < y < yhi
    truth    = (npcoords >= [area[0][0],area[0][1]]) & \
               (npcoords <= [area[1][0],area[1][1]])

    #Check if each pair of values are both in their range together
    truth2   = truth[1:,0] & truth[:-1,1]
    count    = truth2.sum()

    return count

def get_phi_psi_hits( phi_psi ):

    motifs = {
    "alpha_in"  : [ ( ( -pi/2, -5*pi/12 ) , ( -pi/6, -pi/12    )  ) ],
    "alpha_in"  : [ ( ( -pi/2, -5*pi/12 ) , ( -pi/6, -pi/12    )  ) ],
    "alpha_l"   : [ ( ( pi/6, pi/12     ) , ( pi/2, 5*pi/18    )  ) ],
    "alpha_out" : [ ( ( -pi, -pi/2      ) , ( 0.0, pi/6        )  ) ],
    "beta"      : [ ( (-pi, pi/2        ) , (-5*pi/9, pi       )  ) ,
                    ( ( -pi  , -pi      ) , ( -5*pi/9, -5*pi/6 )  ) ,
                    ( ( 5*pi/6, pi/2    ) , ( pi, pi           )  ) ],
    "ppii"      : [ ( ( -5*pi/9, pi/2   ) , ( 0.0, pi          )  ) ,
                    ( ( -5*pi/9, -pi    ) , ( 0.0, -5*pi/6     )  ) ],
    "gamma"     : [ ( ( -pi  , pi/6     ) , ( 0.0, pi/2        )  ) ],
    "gamma_l"   : [ ( ( pi/3, -2*pi/3   ) , ( 2*pi/3, 0.0      )  ) ]
    }

    hits = {}

    for motif_name in motifs:
        hits[motif_name] = 0

        for area in motifs[motif_name]:


            hits[motif_name] += count_in_area( area, phi_psi )

    hits["alpha_out"] -= hits["alpha_in"]

    return hits

def hbond_angle_term( Nx, Ns_COx, Ns_CHx, Cx, Cs_CHx, Cs_NHx ):
    vec_CO_N = Nx - Ns_COx
    vec_CH_N = Nx - Ns_CHx

    vec_N_H = vec_CO_N/linalg.norm(vec_CO_N) + vec_CH_N/linalg.norm(vec_CH_N)
    vec_N_H = vec_N_H/linalg.norm(vec_N_H)

    vec_CH_C = Cx - Cs_CHx
    vec_NH_C = Cx - Cs_NHx

    vec_C_O = vec_CH_C/linalg.norm(vec_CH_C) + vec_NH_C/linalg.norm(vec_NH_C)
    vec_C_O = vec_C_O/linalg.norm(vec_C_O)

    vec_N_C = Cx - Nx
    vec_N_C = vec_N_C/linalg.norm(vec_N_C)

    costhetaN = dot( vec_N_H, vec_N_C )
    costhetaC = dot( vec_C_O, -vec_N_C )

    if costhetaC > 0.0 and costhetaN > 0.0:
        return costhetaC*costhetaC*costhetaN*costhetaN

    return 0.0

def find_PLUM_hbonds( model, same_chain_allowed=True, threshold=0.5 ):

    #settings
    hb_dist_cut=8.0
    sigma=4.11

    #Create list of NH-
    Ns = []
    Cs = []
    pairs_within_cutoff = []
    hbond_pairs = []
    for atom in model.get_atoms():
        if atom.name == 'N':
            Ns.append(atom)
        elif atom.name == 'C':
            Cs.append(atom)

    #Build list of pairs qualifying on:
    #not pro, within dist cutoff, not adjacent or same residue
    for N in Ns:

        #check N doesn't belong to Pro; can't Hbond
        if N.get_parent().get_resname() == '  P':
            continue

        for C in Cs:

            #Within cutoff?
            if linalg.norm( N.get_coord() - C.get_coord() ) > hb_dist_cut:
                continue

            #Not on adjacent or the same residues?
            dist = abs( N.parent.get_id()[1] - C.parent.get_id()[1] )
            on_same_chain = ( N.parent.parent.get_id() == C.parent.parent.get_id() )

            if on_same_chain == True:
                if same_chain_allowed == False or dist < 2:
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
        dist = linalg.norm( N.get_coord() - C.get_coord() )
        dist_term = 5*(sigma/dist)**12 - 6*(sigma/dist)**10

        pre_energy = -dist_term*angle_term

        if pre_energy > threshold:
            hbond_pairs.append([N_resid, C_resid, pre_energy])

    return hbond_pairs

def get_SC_capts( model, threshold, res_separation ):

    SC_atoms = []
    captured = []

    for atom in model.get_atoms():
        if ['N', 'CA', 'C'].count(atom.name) == 0:
            SC_atoms.append(atom)

    for i1, atom1 in enumerate(SC_atoms):
        for atom2 in SC_atoms[i1+1:]:

            same_chain = ( atom1.parent.parent.get_id() == atom2.parent.parent.get_id() )
            dist = abs( atom1.parent.get_id()[1] - atom2.parent.get_id()[1] )

            if same_chain and dist < res_separation:
                continue

            dist = linalg.norm( atom1.get_coord() - atom2.get_coord() )

            if dist < threshold:
                captured.append( (atom1, atom2) )

    return captured

def resID_to_local_resID(resID, res_per_chain):
    local_resID = resID

    while local_resID > 0:
        local_resID -= res_per_chain

    return local_resID+res_per_chain

def hbond_hits_2D(PDB_files,res_per_chain,same_chain_allowed):

    # PLUM hbonds
    models = 0
    pair_hits_ar = zeros([30,30],dtype=int)

    for PDB_file in PDB_files:
        print PDB_file
        for i_model, model in enumerate(bp.PDBParser().get_structure("",PDB_file)):

            models += 1
            hbond_pairs =  find_PLUM_hbonds( model, same_chain_allowed=same_chain_allowed, threshold = 0.2 )
            pair_hits_this = defaultdict(int, {})

            for N_res, C_res, val in hbond_pairs:

                global_pair = frozenset([N_res, C_res])
                if pair_hits_this[global_pair] == 0:
                    pair_hits_this[global_pair] = 1

                    local_resID1 = resID_to_local_resID(N_res, res_per_chain)
                    local_resID2 = resID_to_local_resID(C_res, res_per_chain)
                    pair_hits_ar[local_resID1-1,local_resID2-1] += 1
                    pair_hits_ar[local_resID2-1,local_resID1-1] += 1

    print pair_hits_ar
    print ""
    print pair_hits_ar.astype(float)/(models*8)

if __name__ == '__main__':

    PDB_files=sys.argv[1:]
    hbond_hits_2D(PDB_files,30,False)

