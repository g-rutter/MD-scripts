#!/usr/bin/env python

import chain_analyse_data
import sys
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.gridspec as gridspec

#plt.figure(figsize = (4,4))
#gs1 = gridspec.GridSpec(4, 4)
#gs1.update(wspace=0.025, hspace=0.05) # set the spacing between axes. 

colorcodes=cycle(['b', 'g', 'r', 'c', 'm', 'y'])

chain = sys.argv[1]

N_res = len(chain)
N_res_range = np.arange(N_res)

fig = plt.figure()
gs1 = gridspec.GridSpec(2,3)
gs1.update(wspace=0.025,hspace=0.05)

fig.suptitle(chain)
ax1 = fig.add_subplot(gs1[0], projection='3d')
ax2 = fig.add_subplot(gs1[1], projection='3d')
ax3 = fig.add_subplot(gs1[2], projection='3d')
ax4 = fig.add_subplot(gs1[3])
ax5 = fig.add_subplot(gs1[4])

xpos = np.repeat( N_res_range, N_res ).flatten()
ypos = np.array( [ N_res_range for i in range(N_res) ] ).flatten()
zpos = np.zeros(N_res**2)

dx = dy = 0.3
e1 = np.empty( [N_res,N_res] )
d1 = np.empty( [N_res,N_res] )
h1 = np.empty( [N_res,N_res] )
e2 = np.empty( [2,N_res] )
d2 = np.empty( [N_res,N_res] )

colors = np.empty ( [N_res,N_res],dtype=str )

for index1, residue1 in enumerate(chain):
    thiscolor = [colorcodes.next()]*N_res
    colors[index1] = np.array(thiscolor)

    e2[0,index1] = chain_analyse_data.well_depth[residue1]['NH']
    e2[1,index1] = chain_analyse_data.well_depth[residue1]['CO']

    d2[0,index1] = chain_analyse_data.well_diam[residue1]['NH']
    d2[1,index1] = chain_analyse_data.well_diam[residue1]['CO']

    for index2,residue2 in enumerate(chain):
        e1[index1,index2] = chain_analyse_data.well_depth[residue1][residue2]
        d1[index1,index2] = chain_analyse_data.well_diam[residue1][residue2]
        h1[index1,index2] = chain_analyse_data.sphere_diam[residue1][residue2]

for ax,zdata,title in zip([ax1,ax2,ax3],[e1,d1,h1],['Well depth (SC-SC)','Well diameter (SC-SC)','Sphere diameter (SC-SC)']):
    ax.bar3d(xpos, ypos, zpos, dx, dy, zdata.flatten(), color=colors.flatten(), alpha=0.5)
    ax.set_xticklabels( list(chain) )
    ax.set_yticklabels( list(chain) )
    ax.set_title(title)

for ax,ydata,title in zip([ax4,ax5],[e2,d2],['Well depth (SC-BB)','Well diameter (SC-BB)']):
    ax.bar(N_res_range, ydata[0], dx, color=colorcodes.next())
    ax.bar(N_res_range+dx, ydata[1], dx, color=colorcodes.next())
    ax.set_xticklabels( list(chain) )
    ax.set_title(title)

e1min = e1.min()
e1max = e1.max()
ax1.set_zbound([e1min,e1max])

plt.show()
