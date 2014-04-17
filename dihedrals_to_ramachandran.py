#!/usr/bin/env python
# encoding: utf-8

#Note minus sign in Psis - could be because PDBParser is calcing it oddly
#types of interpolation: bilinear, nearest, bicubic

import math
import matplotlib.pyplot as plt
import numpy as np
from sys import argv

try:
   fn = argv[1]
except IndexError:
   print 'Defaulting to filename = dihedrals_formatted'
   fn = 'dihedrals_formatted'

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

phis_all = []
psis_all = []
weights  = []

with open(fn) as difile:
   for line in difile:
      numbers = line.split(' ')
      psis_all.append( -float(numbers[0]) )
      phis_all.append( float(numbers[1]) )

#heatmap, xedges, yedges = np.histogram2d(psis_all, phis_all, bins=bins, normed=False, range=[[-180,180]]*2, weights=weights)
heatmap, xedges, yedges = np.histogram2d(psis_all, phis_all, bins=bins, normed=True, range=[[-180,180]]*2)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel("Psi (degrees)")
ax.set_xlabel("Phi (degrees)")

imageplot = ax.imshow(heatmap, extent=[-180.0,180.0]*2)
imageplot.set_interpolation('nearest')
#imageplot.set_clim(0, 300)
fig.colorbar(imageplot)
#plt.show()
plt.savefig('foo.png')
