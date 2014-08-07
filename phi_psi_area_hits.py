#!/usr/bin/env python
# encoding: utf-8

#This script takes a pdb file and calculates the portion of hits in area A,
#for each portion of the chain(s) as specified in the settings.

from math import pi, cos
import Bio.PDB as bp
import numpy as np
import argparse
from sys import argv

#################
#  Definitions  #
#################

area_definitons = (
    #bottomleft
    ( ( -np.pi,  000.0 ),   #phi
      ( -np.pi,  000.0 ) ), #psi

    #topleft
    ( ( -np.pi,  000.0 ),   #phi
      (  000.0,  np.pi ) ), #psi

    #bottomright
    ( (  000.0,  np.pi ),   #phi
      ( -np.pi,  000.0 ) ), #psi

    #topright
    ( (  000.0,  np.pi ),   #phi
      (  000.0,  np.pi ) ), #psi

    #alpha_helix
    ( ( -1.92,  -0.66  ),   #phi
      ( -1.40,  -0.26  ) ), #psi
    )

area_names = ( 'bottomleft', 'topleft', 'bottomright', 'topright',
               'alpha_helix' )

#Full n16n and Tiffs n16N regions
tiffs_regions = ( (0,29), (0,7), (8,15), (16,29) )

#Every single peptide bond
every_peptide = ((0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),
                 (10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,17),
                 (17,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,24),
                 (24,25),(25,26),(26,27),(27,28),(28,29)                      )
##############
#  Settings  #
##############

parser = argparse.ArgumentParser(description='Using a PDB trajectory of a chain\
    with Ramachandran angles, calculate the relative population of two areas of\
    the Ramachandran A and B.')
parser.add_argument('PDB_file',  type=str)

args = parser.parse_args()

#My n16N regions (with alpha-helix disrupters as boundaries)
chain_regions = tiffs_regions

counts     = np.zeros( [ len(chain_regions), len(area_definitons) ], dtype = int)
all_counts = np.zeros( [ len(chain_regions) ] )

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

    #Check if each value is in its range
    #i.e. xlo < x < xhi and ylo < y < yhi
    truth    = (npcoords >= [area[0][0],area[1][0]]) & \
               (npcoords <= [area[0][1],area[1][1]])

    #Check if each pair of values are both in their range together
    truth2   = truth[1:,0] & truth[:-1,1]
    count    = truth2.sum()

    return count

######################################
#  Read in phi_psi angles from traj  #
######################################

for i_model, model in enumerate(bp.PDBParser().get_structure("",args.PDB_file)):

    for chain in model:

        poly = bp.Polypeptide.Polypeptide(chain)
        phi_psi = poly.get_phi_psi_list()

        for i_region,region in enumerate(chain_regions):

            phi_psi_region = phi_psi[region[0]:region[1]+1]
            all_counts[i_region] += len(phi_psi_region)-1

            for j_area, area in enumerate(area_definitons):

                counts[i_region][j_area] += count_in_area(area, phi_psi_region)

############
#  Output  #
############

for i_region, region in enumerate(chain_regions):

    print ""
    print ""
    print "Region ("+str(i_region)+"): ", region
    for j_area, area_name in enumerate(area_definitons):
        print "({0}) {1:12}: {2}".format(j_area, area_names[j_area], float(counts[i_region][j_area])/all_counts[i_region])

    print ""
