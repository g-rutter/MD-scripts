#!/usr/bin/env python

import chain_analyse_data
import sys
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
from mpl_toolkits.mplot3d import Axes3D
from prettytable import PrettyTable

import matplotlib.gridspec as gridspec

colorcodes=cycle(['b', 'g', 'r', 'c', 'm', 'y'])
chain = sys.argv[1]

###################
#  Make datasets  #
###################

N_res = len(chain)

well_depth_SCSC = np.empty( [N_res,N_res] )
well_diam_SCSC  = np.empty( [N_res,N_res] )
hard_diam_SCSC  = np.empty( [N_res,N_res] )
well_depth_BBSC = np.empty( [2,N_res] )
well_diam_BBSC  = np.empty( [2,N_res] )

colors = np.empty ( [N_res,N_res],dtype=str )

for index1, residue1 in enumerate(chain):
    thiscolor = [colorcodes.next()]*N_res
    colors[index1] = np.array(thiscolor)

    well_depth_BBSC[0,index1] = chain_analyse_data.well_depth[residue1]['NH']
    well_depth_BBSC[1,index1] = chain_analyse_data.well_depth[residue1]['CO']

    well_diam_BBSC[0,index1] = chain_analyse_data.well_diam[residue1]['NH']
    well_diam_BBSC[1,index1] = chain_analyse_data.well_diam[residue1]['CO']

    for index2,residue2 in enumerate(chain):
        well_depth_SCSC[index1,index2] = chain_analyse_data.well_depth[residue1][residue2]
        well_diam_SCSC[index1,index2] = chain_analyse_data.well_diam[residue1][residue2]
        hard_diam_SCSC[index1,index2] = chain_analyse_data.sphere_diam[residue1][residue2]

datasets_2D = [well_depth_SCSC, well_diam_SCSC, hard_diam_SCSC]
datasets_1D = [well_depth_BBSC, well_diam_BBSC]

########################
#  Print data summary  #
########################

table = PrettyTable(['Category', 'Min', 'Max', 'Ave'])
categories = ['(SC-SC) Well depth', '(SC-SC) Well diam ', '(SC-SC) Sphere diam', '(SC-BB) Well depth', '(SC-BB) Well diam ']

for setname, dataset in zip( categories, datasets_2D + datasets_1D ):
    table.add_row([ setname, dataset.min(), dataset.max(), dataset.mean()] )
print table

###############################
#  Set up figure environment  #
###############################

fig = plt.figure()
gs1 = gridspec.GridSpec(2,3)
gs1.update(wspace=0.10,hspace=0.15)

fig.suptitle(chain)
ax1 = fig.add_subplot(gs1[0], projection='3d')
ax2 = fig.add_subplot(gs1[1], projection='3d')
ax3 = fig.add_subplot(gs1[2], projection='3d')
ax4 = fig.add_subplot(gs1[3])
ax5 = fig.add_subplot(gs1[4])

N_res_range = np.arange(N_res)
xpos = np.repeat( N_res_range, N_res ).flatten()
ypos = np.array( [ N_res_range for i in range(N_res) ] ).flatten()
zpos = np.zeros(N_res**2)

##################
#  Create plots  #
##################

dx = dy = 0.3

for ax,zdata,title in zip([ax1,ax2,ax3],datasets_2D,['Well depth (SC-SC)','Well diameter (SC-SC)','Sphere diameter (SC-SC)']):
    ax.bar3d(xpos, ypos, zpos, dx, dy, zdata.flatten(), color=colors.flatten(), alpha=0.5)
    ax.set_xticklabels( list(chain) )
    ax.set_yticklabels( list(chain) )
    ax.set_title(title)

for ax,ydata,title in zip([ax4,ax5],datasets_1D,['Well depth (SC-BB)','Well diameter (SC-BB)']):
    ax.bar(N_res_range, ydata[0], dx, color=colorcodes.next())
    ax.bar(N_res_range+dx, ydata[1], dx, color=colorcodes.next())
    ax.set_xticklabels( list(chain) )
    ax.set_title(title)

e1min = well_depth_SCSC.min()
e1max = well_depth_SCSC.max()
ax1.set_zbound([e1min,e1max])

plt.show()
