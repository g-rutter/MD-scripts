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
parser.add_argument('-b', '--bins', type=int, nargs='?', default=100,
                    help='Number of bins for Ramachandran in each dimension.')
parser.add_argument('--nocolour', '-n', action='store_true', default=False)
parser.add_argument('PDB_file',  type=str, nargs='+')

args = parser.parse_args()
print args

bins=args.bins
nocolour=args.nocolour
pdb_files=args.PDB_file

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

for pdb_file in pdb_files:
   for model in bp.PDBParser().get_structure("", pdb_file):
      for chain in model:
         poly = bp.Polypeptide.Polypeptide(chain)

         for phi_psi in poly.get_phi_psi_list():
            #print degrees(phi_psi[0]), degrees(phi_psi[1])
            if phi_psi[0] != None and phi_psi[1] != None:
               phis_all.append( degrees(phi_psi[0]))
               psis_all.append( degrees(-phi_psi[1]))

heatmap, xedges, yedges = np.histogram2d(psis_all, phis_all, bins=bins, normed=True, range=[[-180,180]]*2)

#################
#  Plot in plt  #
#################

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel("Psi (degrees)")
ax.set_xlabel("Phi (degrees)")


if nocolour:
    cmap=matplotlib.colors.ListedColormap(['white','black'])
    bounds=[0,np.finfo(type(heatmap.max())).tiny,1]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    imageplot = ax.imshow(heatmap, extent=[-180.0,180.0]*2, cmap=cmap, norm=norm)

else:
    cmap=None
    norm=None
    imageplot = ax.imshow(heatmap, extent=[-180.0,180.0]*2, cmap=cmap, norm=norm)
    fig.colorbar(imageplot)

imageplot.set_interpolation('nearest')
#imageplot.set_clim(2.0, 200)
plt.savefig('Ramachandran.png')
