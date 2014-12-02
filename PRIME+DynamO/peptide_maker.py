#!/usr/bin/python2

from lxml import etree as ET
import sys
import os
import time
import random
import numpy as np
import subprocess
import re

from mylammps import mylammpsclass
import mylammps

##############
#  Settings  #
##############

#Turn on for full output from LAMMPS and DynamO
#including LAMMPS .dcd files.
debug = True

#temp setting. how much beads are shrunk by for close interactions.
scale = 0.85

###############
#  Functions  #
###############

def gauss():
    return str( random.gauss(0.0, 1.0) )

def mkdir_p(path):
    if not os.path.exists(path):
        os.makedirs(path)

def joinStr(list_like):
    return " ".join(str(item) for item in list_like)

#############################################
###  Read in the command line parameters  ###
#############################################

nonglycine_SC    = list('ACDEFHIKLMNPQRSTVWY')
nonglycine_sites = list(nonglycine_SC) + ['NH', 'CH', 'CO']
sites            = list(nonglycine_sites) + ['G']
nonglycine_names = ['Alanine', 'Cysteine', 'Aspartic Acid', 'Glutamic Acid', 'Phenylalanine', 'Histidine', 'Isoleucine', 'Lysine', 'Leucine', 'Methionine', 'Asparagine', 'Proline', 'Glutamine', 'Arginine', 'Serine', 'Threonine', 'Valine', 'Tryptophan', 'Tyrosine' ] + ['Nitrogen+Hydrogen', 'Carbon+Hydrogen', 'Carbon+Oxygen']

try:
    filename = str(sys.argv[3])
except:
    filename = 'Hall_peptide'

psf_fn = filename + '.psf'
xml_fn = filename + '.xml'

try:
    temperature = sys.argv[2]
except:
    temperature = '1.0'

try:
    sys.argv[1] = sys.argv[1].upper()
    if ( sys.argv[1] == 'S1' ):
        sequence = list('PPPWLPYMPPWS')
    elif ( sys.argv[1] == 'N16N' ):
        sequence = list('AYHKKCGRYSYCWIPYDIERDRYDNGDKKC')
    elif ( sys.argv[1][:5] == 'ALPHA' ):
        sequence = list('ACDEFHIKLMNPQRSTVWY')
    else:
        sequence = list(sys.argv[1])
        valid_input = sites + [ str(number) for number in range(0,10) ]
        assert ( [ residue in valid_input for residue in list(sys.argv[1])] )

        matches = re.finditer("[0-9]+", sys.argv[1])

        for match in matches:

            letter_index = match.start()-1
            number = int( match.group() )

            letter = sys.argv[1][letter_index]
            insertion = letter*number
            sequence[letter_index] = insertion

        sequence = filter( lambda x: x.isalpha(), sequence)
        sequence = list(''.join(sequence))

except (ValueError, IndexError) as e:
    print 'Run as ./peptide_maker.py (sequence) [temperature kT = 1.0] [xml_fn = PRIME_peptide].'
    print ''
    raise

date               = time.strftime('%X %x %Z')
box_size_per_res   = 10.0
box_pad            = 5.0
box_size           = 2*box_pad + len(sequence)*box_size_per_res
n_residues         = len(sequence)
n_bb_sites         = 3*n_residues
n_sc_sites         = n_residues - sequence.count('G')
n_sites            = n_bb_sites + n_sc_sites
gap                          = 0.9
dcd_temp_dir       = "tempDCDs"

expanded_sequence = []
sidechain_IDs = []
for letter in sequence:
    expanded_sequence += ['NH','CH','CO']
    if letter != 'G':
        sidechain_IDs.append(len(expanded_sequence))
        expanded_sequence.append(letter)

nonglycine_expanded_sequence = filter(lambda a: a != 'G', expanded_sequence) #'Real' list, e.g. AGA gives NH,CH,CO,NH,CH,CO,NH,CH,CO,A,A

print 'Sequence:' , ''.join(sequence)
print 'File name:' , xml_fn , '\n'

##############################
#  Get Geometry from LAMMPS  #
##############################

res_space = 3.6
prototype_positions = np.array([ [ 0.   ,  0.   ,  1.   ],
                                 [ 1.364,  0.121,  1.504],
                                 [ 2.366,  0.12 ,  0.376],
                                 [ 1.28 ,  2.   ,  0.   ] ]) \
                      + [ -box_size/2 + box_pad, 0.0, 0.0 ]

if debug:
    mkdir_p(dcd_temp_dir)
    lmp = mylammpsclass()
else:
    lmp = mylammpsclass("", ["-screen", "none"])

#set up box
lmp.command( "units real" )
lmp.command( "dimension 3" )
lmp.command( "boundary m m m" )
lmp.command( "atom_style full" )
lmp.command( "atom_modify sort 0 0" )

#fix NVT, atoms and bonds
lmp.command( "read_data box.lmp" )
lmp.command( "fix NVT all nvt temp 300 0.1 100.0" )
lmp.file("LAMMPS_PRIME20.params")

#lennard jones repulsions
#A_NH_sigma = 3.2
#A_CO_sigma = 3.49
#A_NH_cut = ( 2.0**(1.0/6) ) * A_NH_sigma
#A_CO_cut = ( 2.0**(1.0/6) ) * A_CO_sigma

#ljs_cmds = [joinStr(["pair_coeff  1 20 lj/cut/coul/cut 0.6000", A_NH_sigma, A_NH_cut, 10.0]),
           #joinStr(["pair_coeff  1 22 lj/cut/coul/cut 0.6000", A_CO_sigma, A_CO_cut, 10.0]) ]
#lmp.commands(ljs_cmds)

atom_tally = 0
lmp_coords = np.zeros([0,3])

for i_res, res in enumerate(sequence):
    natoms = lmp.get_natoms()

    #Check atoms have been added correctly
    if atom_tally != natoms:
        print "ERROR: Internal atom tally doesn't match LAMMPS' natoms."
        print "Atom tally:", atom_tally
        print "LAMMPS natoms:", natoms
        exit()

    lmp_coords = np.reshape(lmp_coords, [natoms, 3])
    new_coords = mylammps.create_coords(lmp_coords, i_res, sequence, res_space, prototype_positions)

    #Add next residue
    residue_atoms = [ 'NH', 'CH', 'CO']
    if res != 'G':
        residue_atoms += res

    atom_tally += len(residue_atoms)

    #Create atoms for next residue
    for i_atom, atom in enumerate(residue_atoms):
        coordstr = joinStr(new_coords[natoms+i_atom])
        typestr  = str( nonglycine_sites.index(atom) + 1 )
        create_cmd = "create_atoms " + typestr + " single " + coordstr
        lmp.command( create_cmd )

    #Create groups to use for bond/create and addforce fixes
    mylammps.create_groups(lmp, natoms, residue_atoms, sequence, i_res)

    bondcmds = []

    #BB-BB bonds
    #                      ID     group  bond/create Nevery itype jtype Rmin bondtype
    bondcmds.append( "fix  NH_CH  last1  bond/create 1      20     21   10.0 1" )
    bondcmds.append( "fix  CH_CO  last1  bond/create 1      21     22   10.0 2" )
    bondcmds.append( "fix  CO_NH  CO_NH  bond/create 1      22     20   10.0 3" )

    #BB-BB pseudobonds
    #                      ID     group  bond/create Nevery itype jtype Rmin bondtype
    bondcmds.append( "fix  NH_CO  last1  bond/create 1      20    22    10.0 4" )
    bondcmds.append( "fix  CH_NH  CH_NH  bond/create 1      21    20    10.0 5" )
    bondcmds.append( "fix  CO_CH  CO_CH  bond/create 1      22    21    10.0 6" )
    bondcmds.append( "fix  CH_CH  last2  bond/create 1      21    21    10.0 7" )

    bond_fix_names = [ "NH_CH", "CH_CO", "CO_NH", "NH_CO", "CH_NH", "CO_CH", "CH_CH" ]

    #BB-SC bonds and pseudobonds
    if res != 'G':
        sc_type = nonglycine_sites.index(res) + 1
        NH_SC_type = 5 + sc_type*3
        CH_SC_type = 6 + sc_type*3
        CO_SC_type = 7 + sc_type*3
        #                      ID     group  bond/create Nevery itype   jtype           Rmin     bondtype
        bondcmds.append( "fix  NH_SC  last1  bond/create 1      20    "+str(sc_type)+"  10.0 " + str(NH_SC_type))
        bondcmds.append( "fix  CH_SC  last1  bond/create 1      21    "+str(sc_type)+"  10.0 " + str(CH_SC_type))
        bondcmds.append( "fix  CO_SC  last1  bond/create 1      22    "+str(sc_type)+"  10.0 " + str(CO_SC_type))

        bond_fix_names += ["NH_SC", "CH_SC", "CO_SC"]

    #Add force
    forcecmds = []

    pullstrength = 10.0*i_res

    forcecmds.append( "fix  negativepull  firstbbatom addforce -{0:f} 0.0 0.0". format(pullstrength) )
    forcecmds.append( "fix  positivepull  lastbbatom  addforce  {0:f} 0.0 0.0". format(pullstrength) )
    force_fix_names = ["negativepull", "positivepull"]

    lmp.commands(bondcmds+forcecmds)

    if debug:
        lmp.command( "dump dcd all dcd 20 " + dcd_temp_dir + "/LAMMPS_" + str(i_res) + "_test.dcd" )

    #Run
    lmp.command( "run 1" )
    lmp.minimize()

    lmp.command( "run 2000" ) 
    lmp.minimize()

    lmp.unfixes(bond_fix_names+force_fix_names)

    lmp.command( "run 2000" ) 
    lmp.minimize()

    if debug:
        lmp.command( "undump dcd" )

    lmp_coords = np.array(list( lmp.gather_atoms("x", 1,3) ))

lmp.command( "run 20000" ) 
lmp.minimize()

natoms = lmp.get_natoms()

lmp_coords = np.array(list(lmp.gather_atoms( "x", 1, 3)))
lmp_coords = np.reshape(lmp_coords, [natoms, 3])

lmp.close()

#############################################
###               Set up XML              ###
#############################################

HB_strength = 1.0

DynamOconfig = ET.Element    ( 'DynamOconfig', attrib = {'version' : '1.5.0'} )

Simulation   = ET.SubElement ( DynamOconfig, 'Simulation', attrib = {'lastMFT':"-nan"} )
Properties   = ET.SubElement ( DynamOconfig, 'Properties' )
ParticleData = ET.SubElement ( DynamOconfig, 'ParticleData' )

####ParticleData section####
for ID in range( n_sites ):
    ET.SubElement( ParticleData, 'Pt', attrib = {'ID' : str(ID) } )

    ET.SubElement( ParticleData[ID], 'P', attrib = {"x":str( lmp_coords[ID][0] ), "y":str( lmp_coords[ID][1] ), "z":str( lmp_coords[ID][2] )} )
    ET.SubElement( ParticleData[ID], 'V', attrib = {"x":gauss(), "y":gauss(), "z":gauss()} )

####Simulation section####

Scheduler      = ET.SubElement ( Simulation, 'Scheduler',      attrib = {'Type':'NeighbourList'} )
SimulationSize = ET.SubElement ( Simulation, 'SimulationSize', attrib = dict(zip(['x','y','z'], [str(box_size)]*3)) )
BC             = ET.SubElement ( Simulation, 'BC',             attrib = {'Type':'PBC'} )
Genus          = ET.SubElement ( Simulation, 'Genus')
Topology       = ET.SubElement ( Simulation, 'Topology' )
Interactions   = ET.SubElement ( Simulation, 'Interactions')
Locals         = ET.SubElement ( Simulation, 'Locals' )
Globals        = ET.SubElement ( Simulation, 'Globals' )
SystemEvents   = ET.SubElement ( Simulation, 'SystemEvents' )
Dynamics       = ET.SubElement ( Simulation, 'Dynamics', attrib = {'Type':'Newtonian'} )

Sorter = ET.SubElement ( Scheduler, 'Sorter', attrib = {'Type':'BoundedPQMinMax3'} )

#Add each species and its list of IDs to the XML tree
PRIME_species = ET.SubElement( Genus, 'Species', attrib = {'Mass':'PRIMEData', 'Name':'PRIMEGroups', 'Type':'Point'} )
ET.SubElement( PRIME_species, 'IDRange', attrib = {'Type':'All'} )

temp = ET.SubElement( Globals, 'Global', attrib = {'Type':'Cells','Name':'SchedulerNBList','NeighbourhoodRange':'7.400000000000e+00'})
ET.SubElement( temp, 'IDRange', attrib = {'Type':'All'})

#Interactions section
PRIME = ET.SubElement( Interactions, 'Interaction', attrib = {'Type':'PRIME', 'Name':'Backbone', 'Topology':"PRIMEData", 'HBStrength':str(HB_strength), 'Scale':str(scale)} )
ET.SubElement( PRIME, 'IDPairRange', attrib = {'Type':'All'} )

#Topology section
Structure = ET.SubElement ( Topology, 'Structure', attrib = {'Type':'PRIME', 'Name':'PRIMEData'} )
Molecule  = ET.SubElement ( Structure, 'Molecule', attrib = {'StartID':'0', 'Sequence':''.join(sequence)} )

######################
#  Create PSF files  #
######################

psf_atoms_section = ""
psf_bonds_section = ""

i_res=-1
for i_atom, atom in enumerate(expanded_sequence):

    if atom == 'NH':

        atom_label = 'N'

        i_res+=1
        if sequence[i_res-1] == 'G':
            bond_partner = i_atom - 1
        else:
            bond_partner = i_atom - 2

        #Set equal to itself to signal no bond partner:
        if i_res == 0:
            bond_partner = i_atom

    elif atom == 'CH':
        bond_partner = i_atom-1
        atom_label = 'CA'

    elif atom == 'CO':
        bond_partner = i_atom-1
        atom_label = 'C'

    else:
        bond_partner = i_atom-2
        atom_label = atom+'SC'

    psf_atoms_section += "{0: >8d} {1: <4} {2: <4d} {3: <4} {4: <4} {4: <4} {5: >10} {6: >13} {7: >11}\n".format(i_atom+1, str(0), i_res, sequence[i_res], atom_label, "0.000000", "0.0000", "0")

    if bond_partner != i_atom:
        psf_bonds_section += "{0: >8d}{1: >8d}".format(bond_partner+1,i_atom+1)

        if len(psf_bonds_section) - psf_bonds_section.rfind("\n") > 63:
            psf_bonds_section += "\n"

with open(psf_fn, 'w') as psf_file:
    psf_file.write("PSF\n\n\t1 !NTITLE\n REMARKS " + ''.join(sequence) + " STRUCTURE FILE\n REMARKS DATE: " + date + "\n\n")
    psf_file.write("{0: >8d}".format(n_sites) + " !NATOM\n" + psf_atoms_section + "\n")
    psf_file.write("{0: >8d}".format(n_sites-1) + " !NBOND\n" + psf_bonds_section + "\n\n")

#############################################
###              Write file               ###
#############################################

input_file = open(xml_fn, 'w')
input_file.write('<!-- DynamO input file contains the PRIME20 model of the sequence: ' + ''.join(sequence) + '. -->\n')
input_file.write('<!-- Created on ' +date + '. -->\n')
[ input_file.write(ET.tostring(DynamOconfig, pretty_print=True)) ]
input_file.close()

#Add thermostat and rescale via dynamod:
thermostat_command = [ 'dynamod',  '-T', temperature, '-r', temperature, '-o', xml_fn, '-Z', xml_fn ]
print "Running this command:", " ".join(thermostat_command)
if debug:
    print subprocess.Popen(thermostat_command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
else:
    silent_stdout = subprocess.Popen(thermostat_command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()

#Check config is valid with dynamod:
run_command = ['dynamod', xml_fn, "--check"]
print "Running this command:", " ".join(run_command)
if debug:
    print subprocess.Popen(run_command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
else:
    silent_stdout = subprocess.Popen(run_command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()

if debug:
    #Create trajectory file
    traj_command = ['dynamo2xyz', xml_fn]
    print "Running this command:", " ".join(traj_command)
    with open('traj.xyz', 'w') as trajfile:
        xyz = subprocess.check_output(traj_command)
        trajfile.write(xyz)
