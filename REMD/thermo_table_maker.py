#!/usr/bin/env python3
# encoding: utf-8

#This tool unshuffles a set of LAMMPS dump files (tested on dump local)

import collections
import bisect
from sys import argv
import re

############################################################
#  If no args given, just print out how it should be used  #
############################################################

if len(argv) != 5:
    print ("Run as:", argv[0], "[filestring] [first line] [step_column] [CV_column]")
    print ("[filestring]: Prefix to shuffled logfiles, e.g. log.lammps.")
    print ("[first line]: First line of the log files with the CVs to be tabulated on them.")
    print ("[step_column]: Column of the thermo output which gives the timestep number.")
    print ("[CV_column]: Column of the thermo output which gives the CV.")
    print ("Run in a directory containing the files to be unshuffled and the log.lammps corresponding to them.")
    exit()

swaps        = []
swap_steps   = []
sample_steps = []

filestring  = argv[1]
first_line  = int(argv[2])
step_column = int(argv[3])
CV_column  = int(argv[4])

####################################################
#  function to get temperature index of a replica  #
####################################################

def get_temp ( swaps, swap_steps, step ):
    key = bisect.bisect( swap_steps, step ) - 1

    if key == -1:
        key = 0

    return swaps[key]

def point_files_to_line ( fps, line):
    if type(fps) != list:
        fp  = fps
        fps = [ fp ]

    for fp in fps:
        fp.seek(0)
        for i in range(line-1):
            fp.readline()

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

            if len(swaps) == 0:
                pass
            else:
                print ("Log file goes weird, check out this line:")
                print (logfile_line)

n_replicas = len(swaps[0])

print ("Working with", n_replicas , "replicas.")
print ("Read in", len(swap_steps), "swaps.")

##########################################
#  Check in-files exist, open out-files  #
##########################################

print ("Checking files")

in_fns     = [ filestring + str(i) for i in range(n_replicas) ]
out_fn     = "CVs"
in_files   = []

try:
    in_files = [ open(in_fn, 'r') for in_fn in in_fns ]
except IOError:
    print ("Files", in_fns, "don't all exist. Exiting.")
    exit()

out_file = open(out_fn, 'w')

#####################################
#  Work out stats of the datafiles  #
#####################################

print ("Making list of steps that were sampled.")

point_files_to_line( in_files, first_line)

#Make list of steps that samples were taken on = sample_steps
line = in_files[0].readline()
while line:
    try:
        sample_steps.append ( int(line.split()[step_column-1]) )
    except ValueError:
        break

    line = in_files[0].readline()

point_files_to_line( in_files[0], first_line)
#Check all files start with the same timestep
first_steps = [ int(in_files[i].readline().split()[step_column-1]) for i in range(n_replicas) ]
try:
    assert ( len(set(first_steps)) == 1 )
except AssertionError:
    print("Error: Not all first_steps start with the same timestep on the given line. Here is the list first_steps:\n", first_steps)
    exit()

#Check all files end with the same timestep
smallest_last_step = sample_steps[-1]
for file in in_files[1:]:
    while 1:
        line = file.readline()
        try:
            int(line.split()[step_column-1])
        except (ValueError, IndexError):
            smallest_last_step = min( smallest_last_step, int(prev_line.split()[step_column-1]) )
            break
        prev_line = line

if smallest_last_step < sample_steps[-1]:
    print ("Smallest last step is smaller than the final sample_step. Trimming sample_steps...")
    last_step_i = sample_steps.index(smallest_last_step)
    sample_steps = sample_steps[:last_step_i+1]

print ("Output will give potentials for timesteps", sample_steps[0], "to", sample_steps[-1])

##############
#  Unshuffle  #
###############

print ("Unshuffling")

point_files_to_line( in_files, first_line)

line = ['']*n_replicas

for i in range(len(sample_steps)):
    timestep = sample_steps[i]

    temperature_indices = get_temp( swaps, swap_steps, timestep )
    for j_replica in range(n_replicas):
        temperature_index = temperature_indices[j_replica]
        #full_line = in_files[temperature_index].readline()
        full_line = in_files[j_replica].readline()
        try:
            desired_part = full_line.split()[CV_column-1]
        except IndexError:
            print (in_fns[temperature_index], "has no entry at timestep", timestep)
            exit()
        line[temperature_index] = desired_part

    out_file.write( " ".join(line) )
    out_file.write( "\n" )

#################
#  Close files  #
#################

print ("Closing files")

for file in in_files:
    file.close()
