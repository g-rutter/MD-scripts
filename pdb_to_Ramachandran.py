#!/usr/bin/env python
# encoding: utf-8

#Note minus sign in Psis - could be because PDBParser is calcing it oddly
#types of interpolation: bilinear, nearest, bicubic

import Bio.PDB as bp
import math
import matplotlib
import pickle
#Not interactive, can use over SSH
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import numpy as np
import sys
import argparse

##############
#  Settings  #
##############

parser = argparse.ArgumentParser(description='Make Ramachandran out of pdb trajectory')
parser.add_argument('-b', '--bins', type=int, nargs=1, default=200,
                    help='Number of bins for Ramachandran in each dimension.')
parser.add_argument('--nocolour', '-n', action='store_true', default=False)
parser.add_argument('--skip', '-s', type=int, default=1)
parser.add_argument('--notnormed', action='store_true', default=False)
parser.add_argument('PDB_file',  type=str, nargs='*')

parser.add_argument('--textsize', '-t', type=int, nargs='?', default=23,
                    help='Using 23 as standard in journals.')
parser.add_argument('-p', '--pickle', type=str, nargs=1,
                    help='Pickle the heatmap so it won\'t have to be calculated again. Provide the pickle filename.')

parser.add_argument('-u', '--unpickle', type=str, nargs=1,
                    help='Unpickle the heatmap from the provided filename so it doesn\'t have to be calculated now.')

args = parser.parse_args()
print args

bins=args.bins
skip=args.skip
nocolour=args.nocolour
notnormed=args.notnormed
pdb_files=args.PDB_file

#validation
if len(pdb_files) == 0 and args.unpickle == None:
    print "Either provide PDB file(s) or a Python pickle for the heatmap."
    print "Try -h"
    exit()

if args.unpickle != None and args.pickle != None:
    print "Can't both pickle and unpickle. Pick one!"
    exit()

if len(pdb_files) > 0:
    out_file = 'Ramachandran_' + pdb_files[0] + '.eps'
else:
    out_file = 'Ramachandran_' + args.unpickle[0] + '.eps'

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

phis_all = []
psis_all = []

i=0

for pdb_file in pdb_files:
    print pdb_file
    for model in bp.PDBParser().get_structure("", pdb_file):

        if i%skip == 0:
            for chain in model:
                poly = bp.Polypeptide.Polypeptide(chain)

                for phi_psi in poly.get_phi_psi_list():
                    if phi_psi[0] != None and phi_psi[1] != None:
                        phis_all.append( degrees(phi_psi[0]))
                        psis_all.append( degrees(-phi_psi[1]))
        i+=1

if args.unpickle == None:
    heatmap, xedges, yedges = np.histogram2d(psis_all, phis_all, bins=bins, normed=False, range=[[-180,180]]*2)
else:
    with open(args.unpickle[0]) as unpickle_file:
        heatmap = pickle.load(unpickle_file)

if args.pickle != None:
    with open(args.pickle[0], 'w') as pickle_file:
        pickle.dump(heatmap, pickle_file)

#print heatmap

if notnormed == False:
    heatmap = heatmap/(heatmap.max())

#################
#  Plot in plt  #
#################

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel(r'$\psi$', fontsize=(args.textsize)+5)
ax.set_xlabel(r'$\phi$', fontsize=(args.textsize)+5)

ax.tick_params(axis='x', labelsize=args.textsize, pad=10)
ax.tick_params(axis='y', labelsize=args.textsize, pad=10)

ax.xaxis.set_ticks( [-180, -90, 0, 90, 180] )
ax.yaxis.set_ticks( [-180, -90, 0, 90, 180] )

if nocolour:
    cmap=matplotlib.colors.ListedColormap(['white','black'])
    bounds=[0,np.finfo(type(heatmap.max())).tiny,1]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    imageplot = ax.imshow(heatmap, extent=[-180.0,180.0]*2, cmap=cmap, norm=norm)

else:
    cmap=None
    cmap=plt.get_cmap('jet')
    #cmap.set_under('w')
    norm=None
    imageplot = ax.imshow(heatmap, extent=[-180.0,180.0]*2, cmap=cmap, norm=norm)
    cbar = fig.colorbar(imageplot)
    cbar.ax.tick_params(labelsize=args.textsize)

imageplot.set_interpolation('nearest')

#Set lower and upper bounds of the legend
clims = imageplot.get_clim()
print clims
#imageplot.set_clim(clims[1]*0.01, clims[1])
plt.savefig(out_file)
