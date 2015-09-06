#!/usr/bin/env python
# encoding: utf-8

import Bio.PDB as bp
import sys
from numpy import pi, array, linalg, dot
from collections import defaultdict

PDB_file = sys.argv[1]

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

if __name__ == '__main__':

    #PLUM hbonds
    alpha_hbonds = 0
    models = 0
    for i_model, model in enumerate(bp.PDBParser().get_structure("",PDB_file)):

        print "-------------------"
        print "Model", i_model
        print "-------------------"
        hbond_pairs =  find_PLUM_hbonds( model, same_chain_allowed=True, threshold = 0.2 )
        models += 1
        for N_res, C_res, val in hbond_pairs:
            print "{0:2d} {1:2d}: {2:.2f}%".format(N_res, C_res, 100*val)
            if N_res - C_res == 4:
                alpha_hbonds += 1

    print float(alpha_hbonds)/models, "a-hbs per frame"

    ##phi_psi area hits
    #sum_area_hits = defaultdict(lambda: 0)
    #total_phipsi_pairs = 0
    #for model in bp.PDBParser().get_structure("",PDB_file):

        #for chain in model:

            #phi_psi = bp.Polypeptide.Polypeptide(chain).get_phi_psi_list()
            #total_phipsi_pairs += len(phi_psi)-2

            #area_hits = get_phi_psi_hits( phi_psi )

            #for area in area_hits:
                #sum_area_hits[area] += area_hits[area]

    #sum_area_hits['other'] = total_phipsi_pairs-sum(sum_area_hits.values())
    #print sum_area_hits
