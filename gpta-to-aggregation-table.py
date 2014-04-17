#!/usr/bin/env python
# encoding: utf-8

from sys import argv
from numpy import zeros
from argparse import ArgumentParser
from itertools import combinations

#Set up argument parser
parser = ArgumentParser()
parser.add_argument("filestring", help="Prefix to unshuffled colvar files.")
parser.add_argument("columns", type=int, nargs=2, help="First and last column of minimum distances in colvar files; 1-based.")
parser.add_argument('--max', type=float, nargs='?', help="Maximum value of minimum distance to count as aggregated.", default=5.5)

#Interpret input
args  = parser.parse_args()

filestring = args.filestring
cols = args.columns
threshold = args.max

#######################################################################
#               Determine other preliminary properties                #
#######################################################################

N_columns = cols[1]-cols[0]+1

N_peptides = 0
while 1:
   N_peptides += 1
   N_columns_result = N_peptides*(N_peptides-1)/2
   if (N_columns == N_columns_result):
      break
   if (N_columns < N_columns_result):
      print "Error: Wrong column indices."
      exit()

cluster_possibilities = {
         # 1  2
   2 : [ [ 2, 0 ],  #State 0
         [ 0, 1 ] ] #State 1
,
         # 1  2  3
   3 : [ [ 3, 0, 0],  #State 0
         [ 1, 1, 0],  #State 1
         [ 0, 0, 1] ] #State 2
,
         # 1  2  3  4  5  6
   6 : [ [ 6, 0, 0, 0, 0, 0],  #State 0
         [ 4, 1, 0, 0, 0, 0],  #State 1
         [ 2, 2, 0, 0, 0, 0],  #State 2
         [ 3, 0, 1, 0, 0, 0],  #State 3
         [ 1, 1, 1, 0, 0, 0],  #State 4
         [ 0, 3, 0, 0, 0, 0],  #State 5
         [ 2, 0, 0, 1, 0, 0],  #State 6
         [ 0, 1, 0, 1, 0, 0],  #State 7
         [ 0, 0, 2, 0, 0, 0],  #State 8
         [ 1, 0, 0, 0, 1, 0],  #State 9
         [ 0, 0, 0, 0, 0, 1] ] #State 10
}

cluster_list = cluster_possibilities[N_peptides]
pairs        = list( combinations( range(N_peptides), 2) )
N_states     = len(cluster_list)

#######################################################################
#                            Define funcs                             #
#######################################################################

def get_state( entry, cluster_list, N_peptides ):
   'Return the state that an entry corresponds to'
   groups = [ set([ i ]) for i in range(N_peptides) ]

   #Make a list of direct bonds first
   for i, bonded in enumerate( entry ):
      if bonded:
         groups[pairs[i][0]].add(int(pairs[i][1]))
         groups[pairs[i][1]].add(int(pairs[i][0]))

   #Group as clusters
   groups_old = [ set([ i ]) for i in range(N_peptides) ]
   while groups != groups_old:
      groups_old = list(groups)
      for i in range(N_peptides):
         for j in range(i):
            if len( groups[i].intersection( groups[j] ) ) > 0:
               groups[j] = groups[j].union( groups[i] )
               groups[i] = set()

   #Make a list like: [ N_monomers, N_dimers, N_trimers, ... ]
   counts = [ len(groups[i]) for i in range(N_peptides) ]
   state_list = [ counts.count(i) for i in range(1, N_peptides+1) ]

   state = cluster_list.index(state_list)
   return state

#######################################################################
#                  Read all entries and output files                  #
#######################################################################

T=0
input_files = []
while 1:
   T += 1
   try:
      input_files.append ( open(filestring+str(T), 'r') )
   except IOError:
      N_input_files = T
      break

print "N_input_files:", N_input_files
print "N_peptides:", N_peptides

input_files   = [ open(filestring+str(T), 'r') for T in range(N_input_files)      ]
current_lines = [ input_files[i].readline() for i in range(N_input_files) ] #get past comment line
output_files  = [ open("aggregation_state_"+str(n), 'w') for n in range(N_states) ]

while 1:

   current_lines = [ input_files[i].readline() for i in range(N_input_files) ] #get past comment line
   if current_lines[0] == '':
      break
   bonds         = [ [ float(number)<threshold for number in line.split()[cols[0]-1:cols[1]] ] for line in current_lines ]
   states        = [ get_state( bonds[i], cluster_list, N_peptides ) for i in range (N_input_files) ]

   for T in range(N_input_files):

      for n in range(N_states):
         if states[T] == n:
            output_files[n].write("1 ")
         else:
            output_files[n].write("0 ")

   for n in range(N_states):
      output_files[n].write("\n")

for file in input_files+output_files:
   file.close()
