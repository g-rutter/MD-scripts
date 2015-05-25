#!/usr/bin/env python
# encoding: utf-8

import bz2
import pickle
import glob
import re
import sys
import numpy as np
import math
from numpy.linalg import norm
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from collections import defaultdict
from progressbar import ProgressBar

##############
#  Settings  #
##############

#captures won't be counted unless abs(res1-res2) > min_res_dist
min_res_dist = 3
#captures will be printed if they appear in more than this fraction of maximum poss
thresh = 0.1

########################
#  Get and sort files  #
########################

files=sys.argv[1:]
N_files = len(files)

###########################
#  Get system properties  #
###########################

Chain_seq_pattern = re.compile('<Molecule.*Sequence="([A-Z]+)"')
contents = bz2.BZ2File(files[0]).read()
Chain_seq_match = Chain_seq_pattern.findall(contents)

N_residues = len(Chain_seq_match[0])
N_chains = len(Chain_seq_match)

chain_seq = Chain_seq_match[0]
atoms_per_chain = len(chain_seq)*4 - chain_seq.count('G')

print "N_files:", N_files
print "N_residues:", N_residues
print "N_chains:", N_chains
print "Chain seq:", chain_seq
print "Atoms per chain:", atoms_per_chain

#DEFINITIONS
#pIDs are 0-based
#resIDs are 0-based and SKIP 3 between chains
#chainIDs are 0-based
#local pIDs and local resIDs are as if this were chain 0

def resID_to_chainID(resID):
    period = N_residues+3
    chainID = resID/period
    return chainID

def pID_to_resID(pID):
    local_pID = pID_to_local_pID(pID)
    local_resID = local_pID_to_local_resID(local_pID)

    chainID = pID/atoms_per_chain
    resID = (N_residues+3)*chainID

    resID += local_resID
    return resID

def local_pID_to_local_resID(local_pID):
    local_resID = 0

    atom_counter = 4 - chain_seq[:1].count('G')

    while local_pID >= atom_counter:
        local_resID += 1
        atom_counter = len(chain_seq[:local_resID+1])*4 - chain_seq[:local_resID+1].count('G')

    return local_resID

def pID_to_resname(pID):
    local_pID = pID_to_local_pID(pID)
    local_resID = local_pID_to_local_resID(local_pID)
    resname = local_resID_to_resname(local_resID)
    return resname

def pID_to_local_pID(pID):
    resID = 0
    while pID >= atoms_per_chain:
        resID += N_residues+3
        pID -= atoms_per_chain
    return pID

def pID_to_atom_name(pID):
    local_pID = pID_to_local_pID(pID)
    local_resID = local_pID_to_local_resID(local_pID)

    atom_counter = len(chain_seq[:local_resID])*4 - chain_seq[:local_resID].count('G')
    superlocal_pID = local_pID - atom_counter

    if superlocal_pID < 3:
        return ['NH', 'CH', 'CO'][superlocal_pID]
    elif superlocal_pID == 3:
        return chain_seq[local_resID]

def resID_to_local_resID(resID):
    chainID = resID_to_chainID(resID)
    local_resID = resID - chainID*(N_residues+3)
    return local_resID

def local_resID_to_resname(local_resID):
    return chain_seq[local_resID]

def resID_to_resname(resID):
    local_resID = resID_to_local_resID(resID)
    return local_resID_to_resname(local_resID)

###################
#  Process files  #
###################

print "Processing", N_files, "system snaphots..."

pair_pattern = re.compile('Pair ID1="(\d*)" ID2="(\d*)"')
HBond_pair_pattern = re.compile('NH="(\d*)" CO="(\d*)"')

captures_by_PID = defaultdict(int)
captures_by_resID = defaultdict(int)

interchain_capts_resID = defaultdict(int)
interchain_hbs_resID = defaultdict(int)

progress = ProgressBar(maxval=N_files).start()

for i_file, filename in enumerate(files):

    contents = bz2.BZ2File(filename).read()
    capt_entries = pair_pattern.finditer(contents)

    captures_by_resID_temp = defaultdict(int)
    interchain_capts_resID_temp = defaultdict(int)
    interchain_hbs_resID_temp = defaultdict(int)

    for capt_entry in capt_entries:

        pID1 = int(capt_entry.group(1))
        pID2 = int(capt_entry.group(2))
        resID1 = pID_to_resID(pID1)
        resID2 = pID_to_resID(pID2)
        p1name = pID_to_atom_name(pID1)
        p2name = pID_to_atom_name(pID2)

        if abs(resID1 - resID2) < min_res_dist:
            continue

        captures_by_PID [ frozenset((pID1, pID2)) ] += 1

        if captures_by_resID_temp [ frozenset((resID1, resID2 )) ] == 0:

            captures_by_resID_temp [ frozenset((resID1, resID2 )) ] = 1
            captures_by_resID [ frozenset((resID1, resID2 )) ] += 1

        #interchain interaction counter
        if resID_to_chainID(resID1) != resID_to_chainID(resID2):
            if interchain_capts_resID_temp [ resID1 ] == 0:

                interchain_capts_resID_temp [ resID1 ] = 1
                interchain_capts_resID [ resID1 ] += 1

            if interchain_capts_resID_temp [ resID2 ] == 0:

                interchain_capts_resID_temp [ resID2 ] = 1
                interchain_capts_resID [ resID2 ] += 1

    #################
    #  HB analysis  #
    #################

    hb_entries = HBond_pair_pattern.finditer(contents)

    for hb_entry in hb_entries:

        res_NH = int(hb_entry.group(1))
        res_CO = int(hb_entry.group(2))

        if resID_to_chainID(res_NH) != resID_to_chainID(res_CO):
            if interchain_hbs_resID_temp [ res_NH ] == 0:

                interchain_hbs_resID_temp [ res_NH ] = 1
                interchain_hbs_resID [ res_NH ] += 1

            if interchain_hbs_resID_temp [ res_CO ] == 0:

                interchain_hbs_resID_temp [ res_CO ] = 1
                interchain_hbs_resID [ res_CO ] += 1

    progress.update(i_file+1)

print ""

############
#  Output  #
############

##pID
#for pID1 in range(atoms_per_chain*N_chains):
    #for pID2 in range(atoms_per_chain*N_chains):

        #key = frozenset((pID1,pID2))
        #count = captures_by_PID[key]

        #resID1 = pID_to_resID(pID1)
        #resID2 = pID_to_resID(pID2)

        #if count > int( N_files*thresh):
            #print "{1:2d} {0:2s} ({2:1s}) - ".format( pID_to_atom_name(pID1), resID1, pID_to_resname(pID1)),
            #print "{1:2d} {0:2s} ({2:1s}) - ".format( pID_to_atom_name(pID2), resID2, pID_to_resname(pID2)),
            #print "{0:5.2f}%".format( 100*count/N_files )

print ""

##resID
#for resID1 in range( (N_residues+3)*N_chains ):
    #for resID2 in range( (N_residues+3)*N_chains ):

        #key = frozenset((resID1,resID2))
        #count = captures_by_resID[key]

        #if count > int(N_files*thresh):
            #print "{0:1s} ({1:2d}) - ".format( resID_to_resname(resID1), resID1),
            #print "{0:1s} ({1:2d}) - ".format( resID_to_resname(resID2), resID2),
            #print "{0:4.2f}".format( float(count)/N_files )

#interchain resID
#print "interchain capts by resid"
#for resID in range( (N_residues+3)*N_chains ):
    #count = interchain_capts_resID[ resID ]

    #print "{0:2d} {1:5.2f}".format(resID, float(count)/N_files)

#print "interchain hbs by resid"
#for resID in range( (N_residues+3)*N_chains ):
    #count = interchain_hbs_resID[ resID ]

    #print "{0:2d} {1:5.2f}".format(resID, float(count)/N_files)

###########################
#  Chain capture heatmap  #
###########################

heatmap = np.zeros( [N_residues, N_residues], dtype=float )
heatmap2 = np.zeros( [N_residues, N_residues], dtype=float )

for i_res_global in range(N_residues*N_chains):
    for j_res_global in range(N_residues*N_chains):

        chain1 = resID_to_chainID(i_res_global)
        chain2 = resID_to_chainID(j_res_global)

        if chain1 != chain2:
            capts = captures_by_resID[ frozenset((i_res_global, j_res_global)) ]

            if capts > 0:
                i_res_local = resID_to_local_resID(i_res_global)
                j_res_local = resID_to_local_resID(j_res_global)
                heatmap[i_res_local, j_res_local] += capts

heatmap = heatmap / ((N_chains-1)*N_chains/2 * N_files)

with open('heatmap_pickle', 'w') as heatmap_pickle:
    pickle.dump( heatmap, heatmap_pickle )

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel("Residue 1", fontsize=20)
ax.set_xlabel("Residue 2", fontsize=20)

out_file = "heatmap"

#ax.xaxis.set_ticks( [-180, -120, -60, 0, 60, 120, 180] )
#ax.yaxis.set_ticks( [-180, -120, -60, 0, 60, 120, 180] )

cmap=None
cmap=plt.get_cmap('jet')
#cmap.set_under('w')
norm=None
imageplot = ax.imshow(heatmap, cmap=cmap, norm=norm)
fig.colorbar(imageplot)

imageplot.set_interpolation('nearest')

#Set lower and upper bounds of the legend
clims = imageplot.get_clim()
print clims
#imageplot.set_clim(clims[1]*0.01, clims[1])
plt.savefig(out_file)
