#!/usr/bin/env python
# encoding: utf-8

import bz2
import glob
import re

def getNumber(filename):
    number=int(filename[9:-8])
    return number

files = glob.glob('Snapshot.*.xml.bz2')
files.sort(key=getNumber)
N_snaps = len(files)

HBond_list_pattern = re.compile('HBond.*HBonds', re.S)
HBond_pair_pattern = re.compile('NH="(\d*)" CO="(\d*)"')

HBs = 0
alpha_HBs = 0
antialpha_HBs = 0

for filename in files:

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

avg_HBs           = float(HBs)/N_snaps
avg_alpha_HBs     = float(alpha_HBs)/N_snaps
avg_antialpha_HBs = float(antialpha_HBs)/N_snaps

print avg_HBs, avg_alpha_HBs, avg_antialpha_HBs
