#!/usr/bin/env python3
# encoding: utf-8

#This tool takes a dumpfile of info about hbonds and writes a far reduced file
#which reports at each dump-step whether each pair of peptides are bonded or not
#based on how deep in each others' hbond wells they are.

#correct value is -1.8Kcal/mol (lammps real)

import collections
import bisect
import sys
import re
import time


if len(sys.argv) < 5:
   print ("Run as", sys.argv[0], "[input file] [threshold] [N = atoms per peptide] [M = peptides]")
   exit()

in_file = open(sys.argv[1])
out_file = open("reduced_" + sys.argv[1], 'w')
threshold = float(sys.argv[2])
N = int(sys.argv[3])
M = int(sys.argv[4])

if threshold > 0.0:
   exit("Error! Threshold is positive.")

M_entries = str(int(M*(M-1)/2))

line = ''
EOF = False

while EOF == False:

   ##############################
   #  Extract a single dataset  #
   ##############################

   dataset = ''

   if line == '':
      line = in_file.readline()

   while True:
      dataset += line
      line = in_file.readline()

      if line == 'ITEM: TIMESTEP\n':
         break

      if line == '':
         EOF = True
         break

   ########################################
   #  Create 2D list of bonded molecules  #
   ########################################

   bonds = [ [0 for i in range(M)] for j in range(M) ]
   for dataset_line in dataset.split('\n')[5:-1]:
      linelist = dataset_line.split(' ')
      energy = float(linelist[2])

      if (energy < threshold):
         atom1 = float(linelist[0])-1
         atom2 = float(linelist[1])-1
         mol1 = int(atom1/N)
         mol2 = int(atom2/N)
         bonds[mol1][mol2] = 1
         bonds[mol2][mol1] = 1

   #############################
   #  Write reduced dump file  #
   #############################

   out_file.write('\n'.join(dataset.split('\n')[:3]))
   out_file.write('\n')
   out_file.write(M_entries)
   out_file.write('\n')
   out_file.write("ITEM: ENTRIES mol1 mol2 bonded")
   out_file.write('\n')
   for i in range(M):
      for j in range(i):
         outline = str(i) + " " + str(j) + " " + str(bonds[i][j]) + '\n'
         out_file.write(outline)

