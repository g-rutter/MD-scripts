#!/usr/bin/env python
#encoding: utf-8

import bisect
from argparse import ArgumentParser

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

inputfile_str = 'log.lammps'
out_stats_str = 'REMD_stats.txt'
out_xmgra_str = 'REMD_hist.agr'

inputfile = open(inputfile_str, "r")
out_stats = open(out_stats_str, "w")
out_xmgra = open(out_xmgra_str, "w")

######################
#  Read in RE swaps  #
######################

print "Reading in RE swaps"

n_swaps = file_n_lines(inputfile_str)-3

for i_line, line in enumerate(inputfile):

   i_step = i_line - 3
   line_list = line.split(' ')

   try:
      assert line[0].isdigit()

      swap_steps[i_step] = int(line_list[0])
      swaps[i_step][:]   = line_list[1:]

   except AssertionError:
      if line_list[:2] == ['Running', 'on']:
         n_replicas = int(line_list[2])

         swaps      = np.empty([n_swaps,n_replicas], dtype = int)
         swap_steps = np.empty([n_swaps],            dtype = int)

inputfile.close()

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

#stdev of counts at each temperature
stddevs = counts.std(axis = 1)
normed_stddevs = stddevs/n_swaps

#Maxes for graphing
maxcount    = max(stddevs.max(), counts.max())
maxlencount = len(str(maxcount))
maxleni     = len(str(n_replicas-1))

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

print >>out_stats, ""

#STD
print >>out_stats, '{0:{len}}'.format("STD:", len=maxleni+9),
for temp in range(n_replicas):
   print >>out_stats, '{0:{len}}'.format(int(np.round(stddevs[temp])), len=maxlencount),
print >>out_stats, ""
#Normed STD
print >>out_stats, '{0:{len}}'.format("Normed STD:", len=maxleni+9),
for temp in range(n_replicas):
   print >>out_stats, '{0:{len}}'.format(normed_stddevs[temp], len=maxlencount),
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
#print >>out_xmgra, "@    xaxis  tick major" , str(roundup(maxcount/10,1))
#print >>out_xmgra, "@    yaxis  tick minor ticks 1"
print >>out_xmgra, "@    xaxis  tick major 1"
print >>out_xmgra, "@    xaxis  tick minor ticks 0"
print >>out_xmgra, "@    xaxis  tick major size 1.000000"

for i in range(n_replicas):
   print >>out_xmgra, "@target G0.S"+str(i)
   print >>out_xmgra, "@type xy"
   for j in range(n_replicas):
      print >>out_xmgra, j, counts[j][i]
   print >>out_xmgra, "&"

out_xmgra.close()
