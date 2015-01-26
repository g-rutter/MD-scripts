#!/usr/bin/env python
# encoding: utf-8

#This tool unshuffles a set of LAMMPS DCD files

import bisect
from sys import argv, stdout
import re
from subprocess import call
from os import mkdir, rmdir, remove, devnull, rename, listdir
import numpy

##############
#  Settings  #
##############

LOGFILE = 'log.lammps' # File to read Replex swaps from
TMP_DIR = 'temp'       # Name of tmp dir
DIRLIM  = 10000        # Maximum files to fill tmp directory with at once.
ARGLIM  = 1000         # Max DCD files to ask catdcd to stitch together in one command
                       # this is to avoid "Argument list too long" errors and may constrain memory requirements.

############################################################
#  If no args given, just print out how it should be used  #
############################################################

if len(argv) < 5:
    print "Run as:", argv[0], "[filestring] [F] [N] [T] {B}"
    print "[filestring] If the files are called traj_0.dcd .. traj_N.dcd, filestring is 'traj_'."
    print "[F] First dump step."
    print "[N] Frequency of dump steps."
    print "[T] Which temperature to demultiplex"
    print "{B} Which DCD frame to start with. Typically set this to 1 more than the number of frames so far.  = 1 by default."
    print "Run in a directory containing the files to be unshuffled and the log.lammps corresponding to them."
    exit()

############
#  Set-up  #
############

n_swaps = 0
line = ''

with open(LOGFILE) as REfile:
    for line in REfile:

        try:
            assert line[0].isdigit()

            n_swaps += 1

        except AssertionError:
            pass

line_list = line.split(' ')
n_replicas = len(line_list)-1

swaps = numpy.empty((n_swaps,n_replicas), dtype=int)
swap_steps = []

#Assign and validate DCD dump timesteps
first_dump_step = int(argv[2])
N               = int(argv[3])
T_index         = int(argv[4]) #which T to piece together

try:
    begin_frame = int(argv[5])
    print 'Resuming from frame', begin_frame
except IndexError:
    begin_frame = 1

tmp = float(first_dump_step)/N
if float(int(tmp)) != tmp:
    print 'Error: first_dump_step is wrong because it\'s not divisible by N.'
    exit()

################################################################################
#  function to get list temperature indices for replicas 0-N at given timestep #
################################################################################

def get_temp ( swaps, swap_steps, step ):
    key = bisect.bisect( swap_steps, step ) - 1

    if key == -1:
        key = 0

    return swaps[key]

######################
#  Read in RE swaps  #
######################

print "Reading in REMD swaps from log."

n_excluded_lines=0
excluded_lines=[]

with open(LOGFILE) as REfile:
    for i_line, line in enumerate(REfile):

        line_list = line.split(' ')
        try:
            assert line[0].isdigit()

            swap_steps.append( int(line_list[0]) )

            for j_item, item in enumerate(line_list[1:]):
                swaps[i_line][j_item] = int(item)

        except AssertionError:
            n_excluded_lines+=1
            excluded_lines.append( line )

n_replicas = len(line_list)-1
print n_excluded_lines, "were excluded."

###############################################################
#  Check in-files exist, prepare for tmp files and out-files  #
###############################################################

filestring = argv[1]
in_fns     = [ filestring + str(i) + '.dcd' for i in range(n_replicas) ]
fnull      = open(devnull, "w")

for in_fn in in_fns:
    try:
        with open(in_fn, 'rb'):
            pass
    except IOError:
        print "File", in_fn, "doesn't exist. Exiting."
        exit()

###############
#  Unshuffle  #
###############

try:
    mkdir(TMP_DIR)
except OSError: #Assume this means dir exists. Not foolproof
    pass

print "Reordering frames from", n_replicas, "replicas."
stdout.flush()

sample_steps = range(first_dump_step, swap_steps[-1], N)
tmp_DCDs = []
file_prefix = TMP_DIR+"/T"+str(T_index)+"-"

for first in range( begin_frame - 1, len(sample_steps), DIRLIM ):
    #Extract frames, DIRLIM at a time.
    for i, sample_step in enumerate(sample_steps[first:first+DIRLIM]):

        j = i + first
        #Obtain replica index correpsonding to T_index at given sample_step
        pos_array = get_temp( swaps, swap_steps, sample_step )
        T_pos = numpy.where( pos_array == T_index )[0][0]

        sj = str(j+1)
        tmp_DCDs.append(file_prefix+sj+".dcd")

        call( [ "catdcd", "-o", tmp_DCDs[-1], "-first", sj, "-last", sj, in_fns[T_pos] ], stdout=fnull )

    #Recursively stitch together without a massive argument list
    out_fn = file_prefix
    depth = 0

    while len(tmp_DCDs) > 1:
        out_fn += "i" #Number of 'i's will be 'depth' of recursion
        tmp_DCDs2 = [] #The outputted DCDs which will be combined in the next round

        for i in range (0, len(tmp_DCDs), ARGLIM):

            si = str(i*(ARGLIM**depth)+begin_frame) #Ensures file will be named after first frame it represents
            tmp_DCDs2.append( out_fn+si+".dcd" )
            inext = min( i+ARGLIM, len(tmp_DCDs) )

            call( ["catdcd", "-o", tmp_DCDs2[-1] ] + tmp_DCDs[i:inext], stdout=fnull  )
            for file in tmp_DCDs[i:inext]: #remove files ASAP to clear space
                remove(file)

        tmp_DCDs = list(tmp_DCDs2) #Get ready to combine the new DCDs next time round
        depth += 1

    rename(tmp_DCDs[0], file_prefix+str(begin_frame)+".dcd")
    tmp_DCDs = [file_prefix+str(begin_frame)+".dcd"]

rename(tmp_DCDs[0], filestring + "T" + str(T_index) + "-" + str(begin_frame) + ".dcd")

print "Finished temperature index", T_index

fnull.close()

try:
    rmdir(TMP_DIR)
except OSError:
    pass
