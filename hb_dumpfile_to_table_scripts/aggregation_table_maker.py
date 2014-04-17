#!/usr/bin/env python
# encoding: utf-8

from sys import argv
from numpy import zeros

if len(argv) < 3:
   print "Run as:", argv[0], "[filestring] [N_input_files]"
   print "[filestring]: prefix to reduced, unshuffled hb file, e.g. reduced_hb_T"
   print "[N_input_files]: number of input files"
   exit()

filestring = argv[1]
N_input_files = int(argv[2])
EXTRA_LINES = 5

#######################################################################
#        Determine number of peptides and define cluster dict         #
#######################################################################

cluster_possibilities = {}

                             # 1  2
cluster_possibilities[2] = [ [ 2, 0 ],  #State 0
                             [ 0, 1 ] ] #State 1

                             # 1  2  3
cluster_possibilities[3] = [ [ 3, 0, 0],  #State 0
                             [ 1, 1, 0],  #State 1
                             [ 0, 0, 1] ] #State 2

                             # 1  2  3  4  5  6
cluster_possibilities[6] = [ [ 6, 0, 0, 0, 0, 0],  #State 0
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

#Reading file 0 to determine setup
with open( filestring + '0', 'r') as f_T0:

   reading_entries = False
   N_entries = 0
   numbers = []

   while 1:
      line = f_T0.readline()

      if reading_entries:
         try:
            numbers += [ int(number) for number in line.split()[:-1] ]
            N_entries += 1
         except ValueError:
            break

      if line[:13] == 'ITEM: ENTRIES':
         reading_entries = True

N_peptides = max(numbers)+1
cluster_list = cluster_possibilities[N_peptides]
N_states = len(cluster_list)

#Determine number of lines to get n_timesteps
with open( filestring + '0', 'r') as f_T0:
   N_lines = 0
   while f_T0.readline():
      N_lines += 1

N_timesteps = N_lines/(N_entries+EXTRA_LINES)

print "\nFile information:"
print "N_timesteps:", N_timesteps
print "N_peptides:", N_peptides

#######################################################################
#                            Define funcs                             #
#######################################################################

def get_entry( fp , N_entries):
   'Return the entries of the next timestep in file fp as a string'
   lines_to_read = N_entries + EXTRA_LINES

   lines = []

   for i in range(lines_to_read):
      lines.append( fp.readline().strip() )

   return '\n'.join(lines[5:])

def get_state( entry, cluster_list, N_peptides ):
   'Return the state that an entry corresponds to'
   groups = [ set([ i ]) for i in range(N_peptides) ]

   #Make a list of direct bonds first
   for line in entry.split('\n'):
      numbers = line.split(' ')
      if numbers[-1] == '1':
         groups[int(numbers[0])].add(int(numbers[1]))
         groups[int(numbers[1])].add(int(numbers[0]))

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

input_files  = [ open(filestring+str(T), 'r') for T in range(N_input_files)      ]
output_files = [ open("aggregation_state_"+str(n), 'w') for n in range(N_states) ]

entry = ( get_entry( input_files[T], N_entries) )
#N_timesteps = 3

for i in range(N_timesteps):

   for T in range(N_input_files):
      entry = ( get_entry( input_files[T], N_entries) )
      state = ( get_state( entry, cluster_list, N_peptides ) )

      for n in range(N_states):
         if state == n:
            output_files[n].write("1 ")
         else:
            output_files[n].write("0 ")

   for n in range(N_states):
      output_files[n].write("\n")

for file in input_files+output_files:
   file.close()
