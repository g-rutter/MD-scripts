#!/usr/bin/env python
# encoding: utf-8

import pickle
import sys
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from prettytable import PrettyTable
from collections import defaultdict
from progressbar import ProgressBar

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

#Get data
data_source = sys.argv[1]

try:
    outname = sys.argv[2]
except IndexError:
    outname = "heatmap.eps"


heatmap = np.loadtxt(data_source)
# heatmap = heatmap/(10150*8) #<- n16N-8
# heatmap = heatmap/(9861*8) #<- n16NN-8

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel("Residue 1", fontsize=22, labelpad=35)
ax.set_xlabel("Residue 2", fontsize=22, labelpad=35)
plt.gcf().subplots_adjust(bottom=0.15,left=0.25)

#Axes, labels and label sizes
loc = plt.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
ax.yaxis.set_major_locator(loc)
ax.xaxis.set_major_locator(loc)

ax.xaxis.set_tick_params(width=0,labelsize=13)
ax.yaxis.set_tick_params(width=0,labelsize=11)

ax.get_yaxis().set_major_formatter(plt.FixedFormatter( list('AAYHKKCGRYSYCWIPYDIERDRYDNGDKKC') ))
ax.get_xaxis().set_major_formatter(plt.FixedFormatter( list('AAYHKKCGRYSYCWIPYDIERDRYDNGDKKC') ))

#Colour scheme
cmap=None
base_cmap=plt.get_cmap('jet')
cmap = truncate_colormap(base_cmap, 0.0, 0.92)
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
plt.savefig(outname)
