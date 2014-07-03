#!/usr/bin/env python
# encoding: utf-8

#This script takes a pdb file and calculates the ratio of hits in area A to
#hits in area B, for each portion of the chain(s) as specified in the settings.

from math import pi, cos
import Bio.PDB as bp
import numpy as np
from sys import argv

######################
#  Area definitions  #
######################

bottomleft = ( ( -np.pi,  000.0 ),  #phi
               ( -np.pi,  000.0 ) ) #psi

topleft    = ( ( -np.pi,  000.0 ),  #phi
               (  000.0,  np.pi ) ) #psi

alpha_helix = ( ( -1.92,  -0.66  ),  #phi
                ( -1.40,  -0.26  ) ) #psi

##############
#  Settings  #
##############

try:
    pdb_file = argv[1]
except IndexError:
    print "Run as", argv[0], "[pdb_file]"
    exit()

area_a = bottomleft
area_b = topleft

chain_regions = ( (0,29), )

counts = np.zeros( [ len(chain_regions), 2 ] )

#################
#  Function(s)  #
#################

def degrees(rad_angle) :
    """Converts any angle in radians to degrees.

    If the input is None, then it returns None.
    For numerical input, the output is mapped to [-180,180]
    """

    if rad_angle is None :
        return None

    angle = rad_angle * 180 / pi

    while angle > 180 :
        angle = angle - 360

    while angle < -180 :
        angle = angle + 360

    return angle

def count_in_area(area, coordinates):
    """
    area is an area specified as ( ( x_low, x_high ), ( y_low, y_high ) )
    coordinates is an iterable containing N objects, each of which is an
    iterable containing two items forming a 2D coordinate.

    Returns the number of those coordinates which fall within the area.
    """

    npcoords = np.array(coordinates)

    #Check if each value is in its range e.g. xlo < x < xhi
    truth    = (npcoords >= [area[0][0],area[1][0]]) & \
               (npcoords <= [area[0][1],area[1][1]])

    #Check if each pair of values are both in their range together
    truth2   = truth[:,0] & truth[:,1]
    count    = truth2.sum()

    return count

######################################
#  Read in phi_psi angles from traj  #
######################################

for model in bp.PDBParser().get_structure("",pdb_file):

    for chain in model:

        poly = bp.Polypeptide.Polypeptide(chain)
        phi_psi = poly.get_phi_psi_list()

        for i_region,region in enumerate(chain_regions):

            phi_psi_region = phi_psi[region[0]:region[1]]

            counts_a             = count_in_area(area_a, phi_psi_region)
            counts_b             = count_in_area(area_b, phi_psi_region)

            counts[i_region][0] += counts_a
            counts[i_region][1] += counts_b

print counts
exit()
