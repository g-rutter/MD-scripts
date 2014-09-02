#!/usr/bin/env python
# encoding: utf-8

#Note minus sign in Psis - could be because PDBParser is calcing it oddly
#types of interpolation: bilinear, nearest, bicubic

import Bio.PDB as bp
import math
import matplotlib.pyplot as plt
import numpy as np
import sys

bins=100
#plt.colorbar()

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

#n_residues = len ( bp.Polypeptide.Polypeptide ( bp.PDBParser().get_structure("", sys.argv[1])[0].get_list()[0] ).get_sequence() )

pdb_files = filter ( lambda x: x != 'sep', sys.argv[1:] )
separate  = any ( filter ( lambda x: x == 'sep', sys.argv[1:] ) )

#print n_residues
#print separate

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

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel("Psi (degrees)")
ax.set_xlabel("Phi (degrees)")

imageplot = ax.imshow(heatmap, extent=[-180.0,180.0]*2)
imageplot.set_interpolation('nearest')
#imageplot.set_clim(2.0, 200)
fig.colorbar(imageplot)
plt.show()
