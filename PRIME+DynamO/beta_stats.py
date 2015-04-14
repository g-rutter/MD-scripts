#!/usr/bin/env python
# encoding: utf-8

import bz2
import glob
import re
import sys
import numpy as np
from numpy.linalg import norm
from prettytable import PrettyTable
from progressbar import ProgressBar
from collections import defaultdict

##############
#  Settings  #
##############

parallel_max_rad = 1.0472
antipara_min_rad = parallel_max_rad*2
min_bonding_HBs = 3

########################
#  Get and sort files  #
########################

try:
    start_percent=float(sys.argv[1])
except IndexError:
    start_percent=0

def getNumber(filename):
    number=int(filename[9:-8])
    return number

files = glob.glob('Snapshot.*.xml.bz2')
files= filter( lambda file: 'output' not in file, files)

files.sort(key=getNumber)
start_pos=int(start_percent*len(files)/100.0)
N_snaps = len(files)-start_pos

###########################
#  Get system properties  #
###########################

Chain_seq_pattern = re.compile('<Molecule.*Sequence="([A-Z]+)"')
contents = bz2.BZ2File(files[0]).read()
Chain_seq_match = Chain_seq_pattern.findall(contents)
N_residues = len(Chain_seq_match[0])
N_chains = len(Chain_seq_match)

print "N_residues:", N_residues
print "N_chains:", N_chains

def resID_to_chainID(resID):
    period = N_residues+3
    chainID = resID/period
    return chainID

def getVector(chainID, contents):
    startID = chainID*N_residues*4
    LCa_ID  = startID + 5
    ACa_ID  = startID + 21

    LCa_pattern = re.compile('<Pt ID="'+str(LCa_ID)+'">\n.*x="(.*)" y="(.*)" z="(.*)"/')
    ACa_pattern = re.compile('<Pt ID="'+str(ACa_ID)+'">\n.*x="(.*)" y="(.*)" z="(.*)"/')

    LCa_xlist=LCa_pattern.findall(contents)[0]
    ACa_xlist=ACa_pattern.findall(contents)[0]

    LCa_ar = np.array(LCa_xlist, dtype=float)
    ACa_ar = np.array(ACa_xlist, dtype=float)

    return ACa_ar - LCa_ar

def getAngle(vec1, vec2):
    c = np.dot(vec1, vec2)/(norm(vec1)*norm(vec2)) #cosine
    angle = np.arccos( np.clip(c,-1,1) )
    return angle

###################
#  Process files  #
###################

HBond_list_pattern = re.compile('HBond.*HBonds', re.S)
HBond_pair_pattern = re.compile('NH="(\d*)" CO="(\d*)"')

print "Processing", N_snaps, "system snaphots..."
#pbar = ProgressBar(maxval=N_snaps).start()

for i_file, filename in enumerate(files[start_pos:]):

    print "Processing file", filename

    N_anti, N_para = 0, 0

    contents = bz2.BZ2File(filename).read()
    HBond_list_match = HBond_list_pattern.search(contents)

    if HBond_list_match is not None:

        HBond_list = HBond_list_match.group(0)
        entries = HBond_pair_pattern.finditer(HBond_list)
        bonded_peptides = defaultdict(int)

        for entry in entries:
            N_resID = int( entry.group(1) )
            C_resID = int( entry.group(2) )
            N_chainID = resID_to_chainID(N_resID)
            C_chainID = resID_to_chainID(C_resID)

            if N_chainID != C_chainID:
                bonded_peptides[frozenset((N_chainID, C_chainID))] += 1

        for bonded_chain_pair, N_HBs in bonded_peptides.iteritems():
            bonded_chain_pair_set = set(bonded_chain_pair)
            chainID1 = bonded_chain_pair_set.pop()
            chainID2 = bonded_chain_pair_set.pop()
            if N_HBs > min_bonding_HBs:
                vec1 = getVector(chainID1, contents)
                vec2 = getVector(chainID2, contents)
                angle = getAngle(vec1, vec2)

                if angle < parallel_max_rad:
                    N_para += 1
                elif angle > antipara_min_rad:
                    N_anti += 1

    N_tot = N_para + N_anti
    print "Total beta bonded peptides:", N_tot
    if N_tot != 0:
        print "Para: {0:.2f}%\nAnti: {1:.2f}%".format(100.0*N_para/N_tot, 100.0*N_anti/N_tot)
    #pbar.update(i_file+1)

