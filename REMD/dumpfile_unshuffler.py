#!/usr/bin/env python3
# encoding: utf-8

#This tool unshuffles a set of LAMMPS dump files (tested on dump local)

import collections
import bisect
import sys
import re

############################################################
#  If no args given, just print out how it should be used  #
############################################################

if len(sys.argv) == 1:
   print ("Run as:", sys.argv[0], "[filestring e.g. dihedrals_]")
   print ("Run in a directory containing the files to be unshuffled and the log.lammps corresponding to them.")
   exit()

swaps = []
swap_steps = []

####################################################
#  function to get temperature index of a replica  #
####################################################

def get_temp ( swaps, swap_steps, step ):
   key = bisect.bisect( swap_steps, step ) - 1

   if key == -1:
      key = 0

   return swaps[key]

######################
#  Read in RE swaps  #
######################

with open('log.lammps') as REfile:
   for logfile_line in REfile:

      line_list = logfile_line.split(' ')
      try:
         assert logfile_line[0].isdigit()

         swap_steps.append( int(line_list[0]) )
         swaps.append( [ int(item) for item in line_list[1:] ] )

      except AssertionError:
         if line_list[:2] == ['Running', 'on']:
            n_replicas = int(line_list[2])

print ("Reordering", n_replicas , "replicas.")

##########################################
#  Check in-files exist, open out-files  #
##########################################

print ("Checking files.")

filestring = sys.argv[1]
in_fns     = [ filestring + str(i) for i in range(n_replicas) ]
out_fns    = [ filestring + "T" + str(i) for i in range(n_replicas) ]
in_files   = []

try:
   in_files = [ open(in_fn, 'r') for in_fn in in_fns ]
except IOError:
   print ("Files", in_fns, "don't all exist. Exiting.")
   exit()

out_files = [ open(out_fn, 'w') for out_fn in out_fns ]

#####################################
#  Work out stats of the datafiles  #
#####################################

print ("Making list of steps that were sampled.")

test_file = open(in_fns[0], 'r')
sample_steps = []
nextline = False

#Make list of steps that samples were taken on = sample_steps
test_line = test_file.readline()
while test_line:

   if nextline == True:
      sample_steps.append( int(test_line) )
      nextline = False

   if test_line == 'ITEM: TIMESTEP\n':
      nextline = True

   test_line = test_file.readline()

test_file.close()

###############
#  Unshuffle  #
###############

print ("Unshuffling")

line = ['' for i in range(n_replicas) ]

for i, timestep in enumerate(sample_steps):

   for j in range(n_replicas):
      temperature_index = get_temp( swaps, swap_steps, timestep )[j]

      #Need to get the ith entry from jth file.
      #Need to change this to add to dataset until line == 'ITEM:TIMESTEP\n'
      dataset = ''

      if i == 0:
         line[j] = in_files[j].readline()

      while line[j] != '': #i.e. until EOF
         dataset += line[j]
         line[j] = in_files[j].readline()

         if line[j] == 'ITEM: TIMESTEP\n':
            break

      out_files[temperature_index].write(dataset)

#################
#  Close files  #
#################

print ("Closing files")

for file in out_files+in_files:
   file.close()
