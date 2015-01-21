#!/usr/bin/env python
# encoding: utf-8

import bz2
import glob
import re
from prettytable import PrettyTable
from progressbar import ProgressBar

########################
#  Get and sort files  #
########################

def getNumber(filename):
    number=int(filename[9:-8])
    return number

files = glob.glob('Snapshot.*.xml.bz2')
files.sort(key=getNumber)
N_snaps = len(files)

#######################
#  Get system length  #
#######################

Chain_seq_pattern = re.compile('<Molecule.*Sequence="([A-Z]+)"')
contents = bz2.BZ2File(files[0]).read()
Chain_seq_match = Chain_seq_pattern.search(contents)
N_residues = len(Chain_seq_match.group(1))

###################
#  Process files  #
###################

HBond_list_pattern = re.compile('HBond.*HBonds', re.S)
HBond_pair_pattern = re.compile('NH="(\d*)" CO="(\d*)"')

HBs = 0
alpha_HBs = 0
antialpha_HBs = 0
pi_HBs = 0
antipi_HBs = 0

print "Processing", N_snaps, "system snaphots..."
pbar = ProgressBar(maxval=N_snaps).start()

for i_file, filename in enumerate(files):

    contents = bz2.BZ2File(filename).read()
    HBond_list_match = HBond_list_pattern.search(contents)

    if HBond_list_match is not None:

        HBond_list = HBond_list_match.group(0)
        entries = HBond_pair_pattern.finditer(HBond_list)

        for entry in entries:
            HBs += 1

            dist = int( entry.group(1) ) - int ( entry.group(2) )
            if dist == 4:
                alpha_HBs += 1
            if dist == -4:
                antialpha_HBs += 1

            if dist == 5:
                pi_HBs +=1
            if dist == -5:
                antipi_HBs +=1

    pbar.update(i_file+1)

avg_HBs           = float(HBs)           / ( N_residues * N_snaps)
avg_alpha_HBs     = float(alpha_HBs)     / ( (N_residues-4) * N_snaps)
avg_antialpha_HBs = float(antialpha_HBs) / ( (N_residues-4) * N_snaps)
avg_pi_HBs        = float(pi_HBs)        / ( (N_residues-5) * N_snaps)
avg_antipi_HBs    = float(antipi_HBs)    / ( (N_residues-5) * N_snaps)

############
#  Output  #
############

table = PrettyTable(["HB type" , "Percent"])
table.add_row(['All'           , 100.0*avg_HBs])
table.add_row(['Alpha'         , 100.0*avg_alpha_HBs])
table.add_row(['Antialpha'     , 100.0*avg_antialpha_HBs])
table.add_row(['Pi'            , 100.0*avg_pi_HBs])
table.add_row(['Antipi'        , 100.0*avg_antipi_HBs])
table.float_format['Percent'] = "2.2"

print table
