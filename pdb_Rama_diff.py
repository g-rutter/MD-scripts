#!/usr/bin/env python
# encoding: utf-8

#Note minus sign in Psis - could be because PDBParser is calcing it oddly
#types of interpolation: bilinear, nearest, bicubic

import pickle
import Bio.PDB as bp
import math
import matplotlib
#Not interactive, can use over SSH
matplotlib.use('Agg')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse

##############
#  Settings  #
##############

parser = argparse.ArgumentParser(description='Make Ramachandran out of pdb trajectory')
parser.add_argument('unpickle_file', type=str, nargs=2,
                    help='Unpickle the heatmap from the provided filename so it doesn\'t have to be calculated now.')
parser.add_argument('--textsize', '-t', type=int, nargs='?', default=23,
                    help='Using 23 as standard in journals.')

args = parser.parse_args()
print args

unpickle_files=args.unpickle_file

out_file = 'Rama_diff_' + unpickle_files[0] + 'vs' + unpickle_files[1] + '.eps'

print "Processing difference Ramachandran plot of", unpickle_files[0], "minus", unpickle_files[1]

########################
#  Rads->Degrees func  #
########################

def degrees(rad_angle) :
    """Converts any angle in radians to degrees.

    If the input is None, then it returns None.
    For numerical input, the output is mapped to [-180,180]
    """
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle

#############################################
#  Extract phi,psi from pdbs, make 2D hist  #
#############################################

with open(unpickle_files[0]) as unpickle_file:
    heatmap1 = pickle.load(unpickle_file)
    heatmap1 = heatmap1/(heatmap1.sum())
with open(unpickle_files[1]) as unpickle_file:
    heatmap2 = pickle.load(unpickle_file)
    heatmap2 = heatmap2/(heatmap2.sum())

diff_map = heatmap1 - heatmap2
maximum = max( diff_map.max(), -diff_map.min() )

#################
#  Plot in plt  #
#################

matplotlib.rcParams.update({'font.size': args.textsize})

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel(r'$\psi$', fontsize=(args.textsize)+5)
ax.set_xlabel(r'$\phi$', fontsize=(args.textsize)+5)

ax.tick_params(axis='x', labelsize=args.textsize, pad=10)
ax.tick_params(axis='y', labelsize=args.textsize, pad=10)

ax.xaxis.set_ticks( [-180, -90, 0, 90, 180] )
ax.yaxis.set_ticks( [-180, -90, 0, 90, 180] )

cmap=None
cmap=plt.get_cmap('seismic')
#cmap.set_under('w')
norm=None
imageplot = ax.imshow(diff_map, extent=[-180.0,180.0]*2, cmap=cmap, norm=norm)

cbar = fig.colorbar(imageplot)
cbar.formatter.set_powerlimits((0, 0))
cbar.ax.tick_params(labelsize=args.textsize)

imageplot.set_interpolation('nearest')
imageplot.set_clim(-maximum, maximum)

#Set lower and upper bounds of the legend
clims = imageplot.get_clim()
print clims

plt.gcf().subplots_adjust(bottom=0.15)
plt.savefig(out_file)
