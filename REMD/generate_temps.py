#!/usr/bin/env python
# encoding: utf-8

'''Little script to return a list of temperatures of replicas, based on interpolating a function you give'''
from numpy import linspace, round, zeros, interp
from argparse import ArgumentParser
from sys import argv

parser = ArgumentParser(description='Generate a list of temperatures which'\
                                   +' obey the density function passed in.')
parser.add_argument('N', nargs='?', default=None, type=int,
      help='How many temperatures to generate.')
parser.add_argument('--temps', nargs='+', default=[275.0,350.0],
      help='Temperatures at which the densities given apply.')
parser.add_argument('--densities', nargs='+', default=[1.0, 1.0],
      help='Densities at the specified temperatures.')

if len(argv) == 1:
   parser.parse_args(['-h'])

args = parser.parse_args()

print args
print "N =", args.N
print "Temps =", args.temps
print "Densities =", args.densities
print ""

try:
   assert (len(args.temps) == len(args.densities))
except AssertionError:
   print "ERROR: --temps arg length doesn't match --densities."
   exit()

#settings
samples=100000

linear_temps = linspace(275,350,samples)
density = interp( linear_temps, args.temps, args.densities )
density_norm = (density*(args.N-1))/density.sum()

nonlinear_temps=zeros(args.N)
count=1.0
j=0

for i in range(samples):

   if count >= 1.0:
      nonlinear_temps[j]=linear_temps[i]
      count = 0.0
      j += 1

   count+=density_norm[i]

nonlinear_temps[-1]=linear_temps[-1]
maxlen=0
temp_str_list=[]

print "Temperature list:"
for temp in round(nonlinear_temps, 2):
   temp_str = str(temp)
   print temp_str,

   maxlen = max (len(temp_str), maxlen)
   temp_str_list.append(temp_str)
print ""

print ""
print "Formatted for an input script:"
print "variable t world", " ".join(['{i:<{len}}'.format(i=i, len=maxlen) for i in temp_str_list])
print "variable p world", " ".join(['{i:<{len}}'.format(i=i, len=maxlen) for i in range(args.N)])
