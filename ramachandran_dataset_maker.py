#!/usr/bin/env python
# encoding: utf-8

from re import search
from sys import argv, exit
import os
import numpy as np

######################
#  Define functions  #
######################

def findfloat ( line ):
   try:
      return float(search(' -?\d+\.?\d* $', line).group(0)[1:-1])
   except AttributeError: #Indicates the number is really small and has been expressed as MM.MMMe-EE
      return 0.0

#############################
#  Process input arguments  #
#############################

try:
   infn = argv[1]
   infile = open(infn)
except IndexError:
   print 'Run as', argv[0], '[input filename] {sep}'
   print 'Using the sep keyword will produce a separate Ramachandran dataset for each residue.'
   exit()
except IOError:
   print 'File', infn, 'doesn\'t exist.'
   exit()

###############################
#  Determine file properties  #
###############################

n_timesteps     = 0
n_residues      = 1
lines_per_entry = 0

line = infile.readline()
while line:

   if line == 'ITEM: TIMESTEP\n':
      n_timesteps += 1

   if n_timesteps == 1:
      lines_per_entry += 1

      if line[:2] == '1 ':
         n_residues += 1

   line = infile.readline()

#array to store values; holds all values for one timestep at a time
dihedrals = np.ndarray(shape=(n_residues - 2, 2), dtype=float)

#####################################
#  Determine where to output files  #
#####################################

try:
   assert argv[2] == 'sep'
   outfiles = [ open(infn+"_formatted_"+str(i), 'w') for i in range( n_residues - 2 ) ]
except (AssertionError, IndexError):
   outfile = open(infn+"_formatted", 'w')
   outfiles = [ outfile for i in range( n_residues - 2 ) ]

#########################################
#  Final loop to output processed data  #
#########################################

infile.seek(0)
for i in range(n_timesteps):

   #######################################
   #  Construct single-timestep dataset  #
   #######################################

   dataset = ''
   for k in range(lines_per_entry):
      dataset += infile.readline()

   dataset = dataset[:-1]

   ##############################################
   #  Extract dihedrals into dihedral 2D array  #
   ##############################################

   res_no = -1

   for line in dataset.split('\n'):

      if line[:2] == '1 ': #it's a psi angle

         if res_no >= 0: #remove first one with no matching phi
            dihedrals[res_no][1] = findfloat(line)

         res_no += 1

      if line[:2] == '4 ': #it's a phi angle
         try:
            if res_no < n_residues - 2:
               dihedrals[res_no][0] = findfloat(line)
         except IndexError:
            print res_no

   ###################################
   #  Write dihedrals array to file  #
   ###################################

   for i in range(n_residues - 2):
      outfiles[i].write( str(dihedrals[i][0]) + " " + str(dihedrals[i][1]) + "\n" )
