#!/usr/bin/env python
# encoding: utf-8

#Note minus sign in Psis - could be because PDBParser is calcing it oddly
#types of interpolation: bilinear, nearest, bicubic

import Bio.PDB as bp
import math
import matplotlib
#Not interactive, can use over SSH
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse

##############
#  Settings  #
##############

parser = argparse.ArgumentParser(description='Make Ramachandran out of pdb trajectory')
parser.add_argument('-b', '--bins', type=int, nargs='?', default=200,
                    help='Number of bins for Ramachandran in each dimension.')
parser.add_argument('--skip', '-s', type=int, default=1)
parser.add_argument('PDB_file',  type=str, nargs=2)

args = parser.parse_args()
print args

bins=args.bins
skip=args.skip
pdb_files=args.PDB_file

out_file = 'Rama_diff_' + pdb_files[0] + 'vs' + pdb_files[1] + '.png'

print "Processing difference Ramachandran plot of", pdb_files[0], "minus", pdb_files[1]
print "Bins:", bins
print "Skip:", skip

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

phis_all = [ [], [] ]
psis_all = [ [], [] ]

i=0

for i_pdb, pdb_file in enumerate(pdb_files):
    for model in bp.PDBParser().get_structure("", pdb_file):

        if i%skip == 0:
            for chain in model:
                poly = bp.Polypeptide.Polypeptide(chain)

                for phi_psi in poly.get_phi_psi_list():
                    if phi_psi[0] != None and phi_psi[1] != None:
                        phis_all[i_pdb].append( degrees(phi_psi[0]))
                        psis_all[i_pdb].append( degrees(-phi_psi[1]))
        i+=1


heatmap1, xedges, yedges = np.histogram2d(psis_all[0], phis_all[0], bins=bins, normed=False, range=[[-180,180]]*2)
heatmap2, xedges, yedges = np.histogram2d(psis_all[1], phis_all[1], bins=bins, normed=False, range=[[-180,180]]*2)
#print heatmap

heatmap1 = heatmap1/(heatmap1.sum())
heatmap2 = heatmap2/(heatmap2.sum())

diff_map = heatmap1 - heatmap2

maximum = max( diff_map.max(), -diff_map.min() )

#################
#  Plot in plt  #
#################

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel(r'$\psi$', fontsize=20)
ax.set_xlabel(r'$\phi$', fontsize=20)

ax.xaxis.set_ticks( [-180, -120, -60, 0, 60, 120, 180] )
ax.yaxis.set_ticks( [-180, -120, -60, 0, 60, 120, 180] )

cmap=None
cmap=plt.get_cmap('seismic')
#cmap.set_under('w')
norm=None
imageplot = ax.imshow(diff_map, extent=[-180.0,180.0]*2, cmap=cmap, norm=norm)
cbar = fig.colorbar(imageplot)

cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()

imageplot.set_interpolation('nearest')
imageplot.set_clim(-maximum, maximum)

#Set lower and upper bounds of the legend
clims = imageplot.get_clim()
print clims

plt.savefig(out_file)
