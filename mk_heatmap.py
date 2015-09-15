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

data_source = sys.argv[1]
heatmap = np.loadtxt(data_source)
heatmap = heatmap/(10150*8) #<- n16N-8
# heatmap = heatmap/(9861*8) #<- n16NN-8

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_ylabel("Residue 1", fontsize=15)
ax.set_xlabel("Residue 2", fontsize=15)

out_file = "heatmap.eps"

loc = plt.MultipleLocator(base=1.0) # this locator puts ticks at regular intervals
ax.yaxis.set_major_locator(loc)
ax.xaxis.set_major_locator(loc)

ax.xaxis.set_tick_params(width=0)
ax.yaxis.set_tick_params(width=0)

#ax.get_yaxis().set_major_locator(plt.LinearLocator(15))

ax.get_yaxis().set_major_formatter(plt.FixedFormatter( list('AAYHKKCGRYSYCWIPYDIERDRYDNGDKKC') ))
ax.get_xaxis().set_major_formatter(plt.FixedFormatter( list('AAYHKKCGRYSYCWIPYDIERDRYDNGDKKC') ))
#ax.yaxis.set_ticks( [-180, -120, -60, 0, 60, 120, 180] )

cmap=None
cmap=plt.get_cmap('jet')
#cmap.set_under('w')
norm=None
imageplot = ax.imshow(heatmap, cmap=cmap, norm=norm)
fig.colorbar(imageplot)

imageplot.set_interpolation('nearest')

ax.invert_yaxis()

#Set lower and upper bounds of the legend
clims = imageplot.get_clim()
print clims
imageplot.set_clim(0.0, clims[1])
plt.savefig(out_file)
