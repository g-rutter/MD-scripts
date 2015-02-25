#!/usr/bin/env python
# encoding: utf-8

import bz2
import glob
import re
import sys
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

def getVectorAtoms(chainID, contents):
    startID = chainID*N_residues*4
    LCa_ID  = startID + 5
    ACa_ID  = startID + 21

    LCa_pattern = re.compile('<Pt ID="'+str(LCa_ID)+'">\n +(?:<P.*>\n +)? (<V.*>)')
    ACa_pattern = re.compile('<Pt ID="'+str(ACa_ID)+'">\n +(?:<P.*>\n +)? (<V.*>)')
    print LCa_pattern.findall(contents)
###################
#  Process files  #
###################

HBond_list_pattern = re.compile('HBond.*HBonds', re.S)
HBond_pair_pattern = re.compile('NH="(\d*)" CO="(\d*)"')

print "Processing", N_snaps, "system snaphots..."
pbar = ProgressBar(maxval=N_snaps).start()

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
            print chainID1, chainID2
            if chainID1 != chainID2 and val > 2:
                getVectorAtoms(chainID1, contents)


    pbar.update(i_file+1)

