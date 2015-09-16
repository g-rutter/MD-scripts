#!/usr/bin/env python
# encoding: utf-8

import pickle
import sys
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from collections import defaultdict
from progressbar import ProgressBar

#Get data
data_source = sys.argv[1]
heatmap = np.loadtxt(data_source)
# heatmap = heatmap/(10150*8) #<- n16N-8
heatmap = heatmap/(9861*8) #<- n16NN-8

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel("Residue 1", fontsize=22)
ax.set_xlabel("Residue 2", fontsize=22)

#Axes, labels and label sizes
loc = plt.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
ax.yaxis.set_major_locator(loc)
ax.xaxis.set_major_locator(loc)

ax.xaxis.set_tick_params(width=0,labelsize=15)
ax.yaxis.set_tick_params(width=0,labelsize=14)

ax.get_yaxis().set_major_formatter(plt.FixedFormatter( list('AAYHKKCGRYSYCWIPYDIERDRYDNGDKKC') ))
ax.get_xaxis().set_major_formatter(plt.FixedFormatter( list('AAYHKKCGRYSYCWIPYDIERDRYDNGDKKC') ))

#Colour scheme
cmap=None
cmap=plt.get_cmap('jet')
norm=None

#Make heatmap
imageplot = ax.imshow(heatmap, cmap=cmap, norm=norm)
cbar=fig.colorbar(imageplot)
cbar.ax.tick_params(labelsize=18)
imageplot.set_interpolation('nearest')
ax.invert_yaxis()

#Bounds of legend
clims = imageplot.get_clim()
imageplot.set_clim(0.0, clims[1])

#Save file
out_file = "heatmap.eps"
plt.savefig(out_file)
