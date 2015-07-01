#!/usr/bin/env python
# encoding: utf-8

import Bio.PDB as bp
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Get radius of gyration from all atoms')
parser.add_argument('PDB_file',  type=str, nargs='+')

args = parser.parse_args()
print args

###########################
#  Get system properties  #
###########################

models = bp.PDBParser().get_structure("", args.PDB_file[0])
masses = []

for atom in models[0].get_atoms():
    masses.append ( atom.mass )

n_atoms = len(masses)

masses = np.array(masses)
