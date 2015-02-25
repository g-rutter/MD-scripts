#!/usr/bin/env python
# encoding: utf-8

import bz2
import glob
import re
import sys
from numpy import empty, linalg, arccos, pi
from prettytable import PrettyTable
from progressbar import ProgressBar
from collections import defaultdict

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

print N_residues, N_chains

def resID_to_chainID(resID):
    period = N_residues+3
    chainID = resID/period
    return chainID

def getChainVec(chainID, contents):
    startID = chainID*N_residues*4
    LCa_ID  = startID + 5
    ACa_ID  = startID + 21

    LCa_pattern = re.compile('<Pt ID="'+str(LCa_ID)+'">\n +(?:<P.*>\n +)? <V ([xyz]="-?[0-9\.]+") ?([xyz]="-?[0-9\.]+") ?([xyz]="-?[0-9\.]+")')
    ACa_pattern = re.compile('<Pt ID="'+str(ACa_ID)+'">\n +(?:<P.*>\n +)? <V ([xyz]="-?[0-9\.]+") ?([xyz]="-?[0-9\.]+") ?([xyz]="-?[0-9\.]+")')
    coords_tuple_L = LCa_pattern.findall(contents)[0]
    coords_tuple_A = ACa_pattern.findall(contents)[0]

    coords_A = empty([3], dtype=float)
    coords_L = empty([3], dtype=float)

    for i_axis, axis in enumerate(['x','y','z']):
        for coordstr in coords_tuple_A:
            if axis in coordstr:
               coords_A[i_axis]= float( coordstr[3:-1] )
        for coordstr in coords_tuple_L:
            if axis in coordstr:
               coords_L[i_axis]= float( coordstr[3:-1] )

    vec_L_A = coords_A - coords_L
    return vec_L_A

###################
#  Process files  #
###################

HBond_list_pattern = re.compile('HBond.*HBonds', re.S)
HBond_pair_pattern = re.compile('NH="(\d*)" CO="(\d*)"')

print "Processing", N_snaps, "system snaphots..."
pbar = ProgressBar(maxval=N_snaps).start()

parallel_chains, antiparallel_chains = 0, 0

for i_file, filename in enumerate(files[start_pos:]):

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

            bonded_peptides[frozenset((N_chainID, C_chainID))] += 1

        for key, val in bonded_peptides.iteritems():
            setkey = set(key)
            chainID1 = setkey.pop()
            chainID2 = setkey.pop()
            if chainID1 != chainID2 and val > 2:
                chainVec1 = getChainVec(chainID1, contents)
                chainVec2 = getChainVec(chainID2, contents)

                angle = linalg.linalg.dot(chainVec1, chainVec2)
                angle /= (linalg.norm(chainVec1)*linalg.norm(chainVec2))
                angle = arccos(angle)*180/pi
                if angle < 60:
                    status = "Parallel"
                    parallel_chains += 1
                elif angle > 120:
                    status = "Antiparallel"
                    antiparallel_chains += 1
                else:
                    status = "Unclassified"
                print "{0:2d} {1:2d} {2:s}".format(chainID1, chainID2, status)



    pbar.update(i_file+1)

print "Parallel:", parallel_chains
print "Antiparallel:", antiparallel_chains
