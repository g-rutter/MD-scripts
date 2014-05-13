#!/usr/bin/env python
#encoding: utf-8

import bisect
from argparse import ArgumentParser

from glob import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

swaps = []
swap_steps = []

###############
#  Functions  #
###############

def file_n_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def get_temp ( swaps, swap_steps, step ):
    'Get temperature index of replica'
    key = bisect.bisect( swap_steps, step ) - 1

    if key == -1:
        key = 0

    return swaps[key]

def roundup(x,i):
    'Round integers up so they start with i non-rounded digits for neat graphs'

    x       = float(x)
    digits  = len(str(x))
    roundto = float("1"+"0"*(digits-(i+2)))
    rounded = np.ceil(x/roundto)*roundto
    return int(rounded)

#####################
#  Parse arguments  #
#####################

parser = ArgumentParser(description='Stats and graph about REMD temperature swaps.')
parser.add_argument('--log', nargs='?', default='log.lammps', dest='logfile_str',
        help='Name of lammps log file e.g. log.lammps with history of swaps.')
parser.add_argument('--.in', nargs='?', default=None, dest='infile_str',
        help='The input file from the simulation e.g. n16n-1.in. Set to "auto" to autodetect.')
parser.add_argument('--outstats', nargs='?', default='REMD_stats.txt', dest='out_stats_str',
        help='Filename for outputted stats text file.')
parser.add_argument('--outxmg', nargs='?', default='REMD_hist.agr', dest='out_xmgra_str',
        help='Filename for outputted stats text file.')

args = parser.parse_args()
out_stats = open(args.out_stats_str, "w")
out_xmgra = open(args.out_xmgra_str, "w")

#Look for infile automatically.
if args.infile_str == "auto":
    candidate_log_files = glob('*.in')+glob('../*.in')
    try:
        args.infile_str = candidate_log_files[0]
        print "Using", args.infile_str, "for temperature information"
    except IndexError:
        print "No temperature info will be loaded, please put .in in this dir or parent dir",\
                "if you want to use autodetect."

################################
#  Get properties of log file  #
################################

n_header_lines = 0

with open(args.logfile_str, 'r') as logfile:
    for i_line, line in enumerate(logfile):

        line_list = line.split(' ')

        if not line[0].isdigit():
            n_header_lines += 1

        else:
            n_replicas = len( line.split(' ') ) - 1
            break

n_swaps = file_n_lines(args.logfile_str)-n_header_lines

swaps      = np.empty([n_swaps,n_replicas], dtype = int)
swap_steps = np.empty([n_swaps],            dtype = int)

print ""
print "Stats on log file:"
print "n_replicas: ", n_replicas
print "n_swaps: ", n_swaps
print "n_header_lines: ", n_header_lines
print ""

######################
#  Read in RE swaps  #
######################

print "Reading in RE swaps"

with open(args.logfile_str, 'r') as logfile:
    for i_line, line in enumerate(logfile):

        if i_line < n_header_lines:
            continue

        i_step = i_line - n_header_lines
        line_list = line.split(' ')

        swap_steps[i_step] = int(line_list[0])
        swaps[i_step][:]   = line_list[1:]

#####################
#  Get temperature  #
#####################

if args.infile_str != None:
    Temps=True
    x_label="Temperature (K)"
    Xs = np.empty(n_replicas, dtype=float)
    X_major_ticks=10
    with open(args.infile_str, 'r') as infile:

        for line in infile:
            line_list = line.split()

            if len(line_list) == n_replicas+3 and\
                '0' not in line_list and\
                line_list[0] == 'variable':

                print "Detected temperatures are: ",
                for i, T in enumerate(line_list[3:]):
                    Xs[i] = T
                    print T,
                print ""
else:
    Temps=False
    x_label="Temperature index"
    X_major_ticks=1
    print "Using replica index instead of temperatures."
    Xs = np.empty(n_replicas, dtype=int)
    for i in range(n_replicas):
        Xs[i]=i

###############
#  Unshuffle  #
###############

print "Unshuffling"

#replicas_temps[i][j] temperature index of replica i at step j
replicas_temps = np.empty( [ n_replicas, n_swaps ] )
for temp_i in range(n_replicas):
    replicas_temps[temp_i] = np.nonzero(swaps == temp_i)[1]

##############
#  Graph it  #
##############

#if graph:
    #colors = iter(cm.rainbow(np.linspace(0, 0.9, len(replicas_temps))))
    #for y in replicas_temps:
        #plt.plot(swap_steps, y, color=next(colors), linewidth=1)
    #plt.show()

#######################
#  Average trip time  #
#######################

print "Calculating average trip time."

fulltrips = roundtrips = 0

for replica_temps in replicas_temps:

    false_ar        = np.zeros([n_replicas], dtype = bool)
    all_temps_tally = false_ar.copy()
    maxtemp         = False

    for temp in replica_temps:

        if temp == 0:
            if all_temps_tally.sum() == n_replicas:
                all_temps_tally = false_ar.copy()
                fulltrips +=1

            if maxtemp:
                roundtrips += 1
                maxtemp = False

        all_temps_tally[temp] = True

        if temp+1 == n_replicas:
            maxtemp = True

##############################
#  Histogram and trip times  #
##############################

print "Histogramming"

counts = np.empty([n_replicas,n_replicas], dtype=int)

#histogram (main product: counts arary)
for i, replica_temp in enumerate(replicas_temps):
    for temp in range(n_replicas):
        counts[i][temp] = np.nonzero(replica_temp==temp)[0].shape[0]

try:
    full_trip_time = float(swap_steps[-1])/(float(fulltrips)/n_replicas)
except ZeroDivisionError:
    full_trip_time = float("inf")

try:
    round_trip_time = float(swap_steps[-1])/(float(roundtrips)/n_replicas)
except ZeroDivisionError:
    round_trip_time = float("inf")

########################
#  STD devs and slope  #
########################

#stdev of counts at each temperature
stddevs = counts.std(axis = 1)
normed_stddevs = stddevs/n_swaps
slope = np.empty(n_replicas-2, dtype=float)

for i in range(1, n_replicas-1):
    slope[i-1] = np.absolute( counts[i-1]-counts[i+1]).sum()/(float(Xs[i+1]-Xs[i-1])*n_replicas)

##############################
#  Make out_stats text file  #
##############################

#Maxes for graphing
maxcount    = max(stddevs.max(), counts.max())
maxlencount = len(str(maxcount))
maxleni     = len(str(n_replicas-1))

print >>out_stats, n_replicas , "replicas,", "%.2e steps." % (swap_steps[-1])
print >>out_stats, ""

print >>out_stats, "Total trips made by all replicas:"
print >>out_stats, fulltrips, "\"full trips\", starting at the lowest T, visiting every temperature at least once, returning to the lowest."
print >>out_stats, roundtrips, "\"round trips\", starting at the lowest T, visiting the highest, returning to the lowest."

print >>out_stats, ""
print >>out_stats, "Trip time = number_of_timesteps/(number_of_trips/number_of_replicas)"
print >>out_stats, "Full trip time = %.2e timesteps." % full_trip_time
print >>out_stats, "Round trip time = %.2e timesteps." % round_trip_time

print >>out_stats, ""
print >>out_stats, "Histogram of visits to temperatures:"
print >>out_stats, ""
print >>out_stats, " "*(maxleni + maxlencount + 4), ' '.join([ '{0:>{len}}'.format("T"+str(i), len=maxlencount) for i in range(n_replicas) ])
if Temps:
    print >>out_stats, "Temp (K):  "+" "*(maxleni + maxlencount - 8), ' '.join([ '{0:>{len}}'.format(str(int(Xs[i])), len=maxlencount) for i in range(n_replicas) ])
    print >>out_stats, ""
for i, replica_temp in enumerate(replicas_temps):
    print >>out_stats, "Replica", '{0:{len}}'.format(i, len=maxleni)+":",
    for temp in range(n_replicas):
        print >>out_stats, '{0:{len}}'.format(counts[temp][i], len=maxlencount),
        pass
    print >>out_stats, ""

print >>out_stats, ""
print >>out_stats, "Other stats"

#STD
print >>out_stats, '{0:{len}}'.format("STD:", len=maxleni+9),
for temp in range(n_replicas):
    print >>out_stats, '{0:{len}}'.format(int(np.round(stddevs[temp])), len=maxlencount),

#Slope
print >>out_stats, ""
print >>out_stats, '{0:{len}}'.format("Avg slope:", len=maxleni+9),
print >>out_stats, '{0:>{len}}'.format('-', len=maxlencount),
for temp in range(n_replicas-2):
    print >>out_stats, '{0:{len}.0f}'.format(slope[temp], len=maxlencount),
print >>out_stats, '{0:>{len}}'.format('-', len=maxlencount),
print >>out_stats, ""

##Normed STD
#print >>out_stats, '{0:{len}}'.format("Normed STD:", len=maxleni+9),
#for temp in range(n_replicas):
    #print >>out_stats, '{0:{len}}'.format(normed_stddevs[temp], len=maxlencount),
#print >>out_stats, ""

out_stats.close()

############################
#  Make hist xmgrace file  #
############################

print >>out_xmgra, "# Grace project file"
print >>out_xmgra, "#"

print >>out_xmgra, "@version 50122"
print >>out_xmgra, "@with g0"
print >>out_xmgra, "@    WORLD", str(Xs.min()) + ", 0,", str(Xs.max())+",", str(roundup(maxcount,2))
print >>out_xmgra, "@    TITLE \"Visit count to temperatures by replicas\""
print >>out_xmgra, "@    SUBTITLE \"Each line is one of the", n_replicas, "replicas\""
print >>out_xmgra, "@    XAXIS  label \""+x_label+"\""
print >>out_xmgra, "@    YAXIS  label \"Visits\""
#print >>out_xmgra, "@    xaxis  tick major" , str(roundup(maxcount/10,1))
#print >>out_xmgra, "@    yaxis  tick minor ticks 1"
print >>out_xmgra, "@    XAXIS  tick major", X_major_ticks
print >>out_xmgra, "@    XAXIS  tick minor ticks 0"
print >>out_xmgra, "@    XAXIS  tick major size 1.000000"

for i in range(n_replicas):
    print >>out_xmgra, "@target G0.S"+str(i)
    print >>out_xmgra, "@type xy"
    for j in range(n_replicas):
        print >>out_xmgra, Xs[j], counts[j][i]
    print >>out_xmgra, "&"

if Temps:

    print >>out_xmgra, "@    s"+str(n_replicas)+" symbol 4"
    print >>out_xmgra, "@    s"+str(n_replicas)+" symbol size 0.730000"
    print >>out_xmgra, "@target G0.S"+str(n_replicas)
    print >>out_xmgra, "@type xy"
    for j in range(n_replicas):
        print >>out_xmgra, Xs[j], 0
    print >>out_xmgra, "&"

out_xmgra.close()
