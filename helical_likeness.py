#!/usr/bin/env python
# encoding: utf-8

from math import pi, cos
import Bio.PDB as bp
import numpy as np
import cPickle as pickle
from sys import argv

try:
   pdb_file = argv[1]
   helix_frame = argv[2]
except IndexError:
   print "Run as", argv[0], "[pdb_file] [helix_frame/PLUM/ATOM]"
   exit()

pickle_name = pdb_file+".pickle"

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

######################################
#  Read in phi_psi angles from traj  #
######################################

try:
   phi_psi = pickle.load( open ( pickle_name, "rb" ) )
except IOError:
   phi_psi = []
   for model in bp.PDBParser().get_structure("",pdb_file):
      for chain in model:
         poly = bp.Polypeptide.Polypeptide(chain)
         phi_psi.append( poly.get_phi_psi_list() )
   pickle.dump( phi_psi, open( pickle_name, "wb" ) )

phi_psi_0_deg = [ ( degrees( angle[0] ), degrees( angle[1] ) ) for angle in phi_psi[0] ]
print phi_psi_0_deg
exit()

##############################
#  Set ideal phi_psi angles  #
##############################

try:
   helix_frame = int(helix_frame)
   phi_psi_ideal_helix = np.zeros(2)
   phi_psi_ideal_helix = np.average(phi_psi[helix_frame][1:-1], 0)
except ValueError:
   if helix_frame.upper() == 'PLUM':
      phi_psi_ideal_helix = np.array([-0.80392073, -1.10673603])

   elif helix_frame.upper() == 'ATOM':
      phi_psi_ideal_helix = np.array([-1.3086983, -0.55439015])

   elif helix_frame.upper() == 'ALL':
      phi_psi_ideal_helix = np.zeros(2)
      phi_ideal = psi_ideal = 0
      for helix_frame in range(len(phi_psi)):
         phi_psi_ideal_helix = np.average(phi_psi[helix_frame][1:-1], 0)
         phi_ideal += phi_psi_ideal_helix[0]
         psi_ideal += phi_psi_ideal_helix[1]
      phi_ideal /= len(phi_psi)
      psi_ideal /= len(phi_psi)
      phi_psi_ideal_helix = np.array([phi_ideal, psi_ideal])

   else:
      print "Can't understand helix_frame setting", helix_frame
      exit()

print "Using phi,psi ideal helix of" , phi_psi_ideal_helix

#####################################
#  Frame by frame helix similarity  #
#####################################

helix_similarity = np.zeros(len(phi_psi), np.double)

for i, frame in enumerate(phi_psi):
   for pair in frame[1:-1]:
      helix_similarity[i] += 1+cos( pair[0] - phi_psi_ideal_helix[0] )
      helix_similarity[i] += 1+cos( pair[0] - phi_psi_ideal_helix[0] )

   helix_similarity[i] /= 4
   helix_similarity[i] /= len(frame[1:-1])

print "Last 20 values: ", helix_similarity[-20:]

print "Average over all frames", np.average(helix_similarity)
