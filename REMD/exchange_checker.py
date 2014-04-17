#!/usr/bin/env python
#encoding: utf-8

import collections
import bisect
import sys
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from sys import argv

swaps = []
swap_steps = []

###############
#  Functions  #
###############

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

inputfile_str = 'log.lammps'
out_stats_str = 'REMD_stats.txt'
out_xmgra_str = 'REMD_hist.agr'

inputfile = open(inputfile_str, "r")
out_stats = open(out_stats_str, "w")
out_xmgra = open(out_xmgra_str, "w")

#graph = args.graph

######################
#  Read in RE swaps  #
######################

for line in inputfile:

   line_list = line.split(' ')
   try:
      assert line[0].isdigit()

      swap_steps.append( int(line_list[0]) )
      swaps.append( [ int(item) for item in line_list[1:] ] )

   except AssertionError:
      if line_list[:2] == ['Running', 'on']:
         n_replicas = int(line_list[2])

inputfile.close()

###############
#  Unshuffle  #
###############

replicas_temps = [ [] for i in range(n_replicas) ]

for swap in swaps:
   for i in range(n_replicas):
      replicas_temps[i].append(swap.index(i))

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


fulltrips = roundtrips = 0

for replica_temp in replicas_temps:

   all_temps = [ i for i in range (1, n_replicas) ]
   extreme_temps = [ 0, n_replicas-1 ]

   for temp in replica_temp:

      if all_temps == []:
         if temp == 0:
            all_temps = [ i for i in range (1, n_replicas) ]
            fulltrips += 1

      try:
         all_temps.remove(temp)
      except ValueError:
         pass

      if extreme_temps == []:
         extreme_temps = [ 0, n_replicas-1 ]
         roundtrips += 1

      if temp == extreme_temps[0]:
         extreme_temps.remove(temp)

##############################
#  Histogram and trip times  #
##############################

maxleni = maxlencount = maxcount = 0
counts = [ [] for i in range(n_replicas) ]

for i, replica_temp in enumerate(replicas_temps):
   maxleni = max(len(str(i)), maxleni)
   for temp in range(n_replicas):
      counts[i].append(replica_temp.count(temp))
      maxlencount = max(len(str(counts[i][-1])), maxlencount)
      maxcount = max(counts[i][-1], maxcount)

try:
   full_trip_time = float(swap_steps[-1])/(float(fulltrips)/n_replicas)
except ZeroDivisionError:
   full_trip_time = float("inf")

try:
   round_trip_time = float(swap_steps[-1])/(float(roundtrips)/n_replicas)
except ZeroDivisionError:
   round_trip_time = float("inf")

##############################
#  Make out_stats text file  #
##############################

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
for i, replica_temp in enumerate(replicas_temps):
   print >>out_stats, "Replica", '{0:{len}}'.format(i, len=maxleni)+":",
   for temp in range(n_replicas):
      print >>out_stats, '{0:{len}}'.format(counts[temp][i], len=maxlencount),
      pass
   print >>out_stats, ""

out_stats.close()

############################
#  Make hist xmgrace file  #
############################

print >>out_xmgra, "# Grace project file"
print >>out_xmgra, "#"

print >>out_xmgra, "@version 50122"
print >>out_xmgra, "@with g0"
print >>out_xmgra, "@    world 0, 0,", str(n_replicas-1)+",", str(roundup(maxcount,2))
print >>out_xmgra, "@    title \"Visit count to temperatures by replicas\""
print >>out_xmgra, "@    subtitle \"Each line is one of the", n_replicas, "replicas\""
print >>out_xmgra, "@    xaxis  label \"Temperature index\""
print >>out_xmgra, "@    yaxis  label \"Visits\""
print >>out_xmgra, "@    xaxis  tick major" , str(roundup(maxcount/10,1))
print >>out_xmgra, "@    yaxis  tick minor ticks 1"
print >>out_xmgra, "@    xaxis  tick major 10"
print >>out_xmgra, "@    xaxis  tick minor ticks 1"

for i in range(n_replicas):
   print >>out_xmgra, "@target G0.S"+str(i)
   print >>out_xmgra, "@type xy"
   for j in range(n_replicas):
      print >>out_xmgra, j, counts[j][i]
   print >>out_xmgra, "&"

out_xmgra.close()
