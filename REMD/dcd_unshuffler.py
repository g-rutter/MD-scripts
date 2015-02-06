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
    print "[F] First DCD dump step."
    print "[N] Interval between dcd dump steps."
    print "[T] Which temperature to demultiplex"
    print "{B} Which DCD frame to start with. Typically set this to 1 more than the number of frames so far.  = 1 by default."
    print "Run in a directory containing the files to be unshuffled and the log.lammps corresponding to them."
    exit()

############
#  Set-up  #
############

#Assign and validate DCD dump timesteps
first_dcd_dump_step = int(argv[2])
dcd_dump_dt         = int(argv[3])
T_index             = int(argv[4]) #which T to piece together
T_idx_str           = str(T_index)

try:
    begin_frame = int(argv[5])
    print 'Resuming from frame', begin_frame
except IndexError:
    begin_frame = 1

tmp = float(first_dcd_dump_step)/dcd_dump_dt
if float(int(tmp)) != tmp:
    print 'Error: first_dcd_dump_step is wrong because it\'s not divisible by dcd_dump_dt (aka N).'
    exit()

n_swaps_full = 0
first_swap_step = None
line = ''

with open(LOGFILE) as REfile:
    content = REfile.readlines()
    for line in content:

        try:
            if first_swap_step == None:
                assert line[0].isdigit()
                first_swap_step = int( line.split(' ')[0] )

            n_swaps_full += 1

        except AssertionError:
            pass

    t_swap = int( content[-1].split(' ')[0] ) - int ( content[-2].split(' ')[0] )

n_replicas = len( line.split(' ') ) - 1

swap_steps_full  = numpy.array( [ i_swap*t_swap + first_swap_step for i_swap in range(n_swaps_full) ] )
dcd_sample_steps = numpy.array( range(first_dcd_dump_step, swap_steps_full[n_swaps_full-1], dcd_dump_dt) )

###############
#  Functions  #
###############

def get_swap_steps_i ( swap_steps, step ):
    '''
    Returns the position in the swaps array corresponding to the provided step.
    '''

    key = bisect.bisect( swap_steps, step ) - 1

    if key == -1:
        key = 0

    return key

def read_swaps ( LOGFILE, used_swap_steps, T_idx_str ):
    '''
    Produce a numpy array with the replica index of temperature T
    for each used timestep
    '''
    n_swaps_used = len(used_swap_steps)
    T_positions  = numpy.empty( n_swaps_used, dtype = int)

    n_excluded_lines=0
    i_used_swaps = 0

    with open(LOGFILE) as REfile:
        for line in REfile:

            line_list = line.split()

            try:
                assert int(line_list[0]) == used_swap_steps[i_used_swaps]

                T_positions[i_used_swaps] = line_list[1:].index(T_idx_str)
                i_used_swaps += 1

            except (AssertionError, IndexError):
                n_excluded_lines+=1

    n_replicas = len(line_list)-1

    print n_excluded_lines, "lines were excluded."

    return T_positions

######################
#  Read in RE swaps  #
######################

print "Reading in REMD swaps from log."

used_swap_steps = numpy.array( [ swap_steps_full[get_swap_steps_i( swap_steps_full, sample_step)]
                                                                for sample_step in dcd_sample_steps ] )
swap_steps_full = None

T_positions = read_swaps( LOGFILE, used_swap_steps, str(T_index) )

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

tmp_DCDs = []
file_prefix = TMP_DIR+"/T"+str(T_index)+"-"

for first in range( begin_frame - 1, len(dcd_sample_steps), DIRLIM ):
    #Extract frames, DIRLIM at a time.
    for i, sample_step in enumerate(dcd_sample_steps[first:first+DIRLIM]):

        j = i + first
        #Obtain replica index corresponding to T_index at given sample_step
        key = get_swap_steps_i( used_swap_steps, sample_step )
        T_pos = T_positions[key]

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
