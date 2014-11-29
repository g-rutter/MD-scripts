#!/usr/bin/python2

from lxml import etree as ET
import sys
import os
import time
import random
import numpy as np
import subprocess

from mylammps import mylammpsclass
import mylammps

import PRIME20_unbonded
import PRIME20_bonded
import PRIME20_masses

#Turn on for full output from LAMMPS and DynamO
#including LAMMPS .dcd files.
debug = False

def gauss():
    return str( random.gauss(0.0, 1.0) )

def mkdir_p(path):
    if not os.path.exists(path):
        os.makedirs(path)

def joinStr(string, list_like):
    return string.join(str(item) for item in list_like)

def getLambdaDict(well_diam, diam):
    return dict( (k , str( float(well_diam[k]) / float(diam[k]) )) for k in well_diam.keys() )

def comboDict(dicts = []):
    new_dict = {}
    for d in dicts:
        new_dict = dict( new_dict.items() + d.items() )

    return new_dict

#############################################
###       Define interaction dicts        ###
#############################################

attribset              = {'Elasticity':'1'}

#Dicts imported by category
SC_SC_diam      = dict(PRIME20_unbonded.SC_SC_diam)
SC_SC_well_diam = dict(PRIME20_unbonded.SC_SC_well_diam)
SC_SC_lambda    = getLambdaDict(SC_SC_well_diam, SC_SC_diam)
SC_SC_welldepth = dict(PRIME20_unbonded.SC_SC_well_depth)

BB_SC_unbond_diam      = dict(PRIME20_unbonded.BB_SC_diam)
BB_SC_unbond_well_diam = dict(PRIME20_unbonded.BB_SC_well_diam)
BB_SC_unbond_lambda    = getLambdaDict(BB_SC_unbond_well_diam, BB_SC_unbond_diam)
BB_SC_welldepth        = dict(PRIME20_unbonded.BB_SC_well_depth)

BB_SC_bond_diam      = dict(PRIME20_bonded.BB_SC_diam)
BB_SC_bond_well_diam = dict(PRIME20_bonded.BB_SC_well_diam)
BB_SC_bond_lambda    = getLambdaDict(BB_SC_bond_well_diam, BB_SC_bond_diam)

BB_SC_pseudo_diam      = dict(PRIME20_bonded.BB_SC_pseudo_diam)
BB_SC_pseudo_well_diam = dict(PRIME20_bonded.BB_SC_pseudo_well_diam)
BB_SC_pseudo_lambda    = getLambdaDict(BB_SC_pseudo_well_diam, BB_SC_pseudo_diam)

#Combined dicts
unbond_diam       = comboDict( [SC_SC_diam, BB_SC_unbond_diam] )
unbond_lambda     = comboDict( [SC_SC_lambda, BB_SC_unbond_lambda] )
unbond_welldiam   = comboDict( [SC_SC_well_diam, BB_SC_unbond_well_diam] )
welldepth         = comboDict( [SC_SC_welldepth, BB_SC_welldepth] )

#Close (i.e. bonded via neighs) interactions
close_diam   = dict(PRIME20_unbonded.close_diam)
close_lambda = getLambdaDict(unbond_welldiam, unbond_diam)

#pseudobond_lambda      = '1.0450' #}http://www.sciencedirect.com/science/article/pii/S0022283611013672 Materials and Methods
#bond_lambda            = '1.0450' #}

#Interactions will end up written as XML entries to these dicts
unbondint = {}
closeunbondint = {}
pseudobondint = {}
bondint = {}

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
    sys.argv[1] = [ letter.upper() for letter in sys.argv[1] ]
    if ( sys.argv[1] == list('N16N') ):
        sequence = list('AYHKKCGRYSYCWIPYDIERDRYDNGDKKC')
    elif ( sys.argv[1][:5] == list('ALPHA') ):
        sequence = list('ACDEFHIKLMNPQRSTVWY')
    else:
        sequence = list(sys.argv[1])
        assert [ sites.index(residue) for residue in sequence]
except ValueError:
    sys.exit('Run as ./peptide_maker.py (sequence) [temperature kT = 1.0] [xml_fn = Hall_peptide].')

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

#############################################
###      Geometry of one amino acid       ###
#############################################

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

lmp.command( "units real" )
lmp.command( "dimension 3" )
lmp.command( "boundary m m m" )
lmp.command( "atom_style molecular" )
lmp.command( "atom_modify sort 0 0" )

lmp.command( "read_data box.lmp" )
lmp.file("LAMMPS_PRIME20.params")

lj_sigma = 3.0
lj_distance_cutoff = ( 2.0**(1.0/6) ) * lj_sigma

lmp.command( "pair_style table linear 2" )
lmp.command( "pair_coeff * * LAMMPS_repulsive_pair_table.txt repulsive_linear" )

#lmp.command( "pair_style lj/cut 10.0" )
#lmp.command( "pair_coeff * * 0.0000 " + str(lj_sigma) )
#lmp.command( "pair_modify shift yes" )

lmp.command( "special_bonds lj 0 1 1" )

lmp.command( "fix NVT all nvt temp 300 0.1 100.0" )

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
        coordstr = joinStr(" ", new_coords[natoms+i_atom])
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
###             Set up classes            ###
#############################################

class Species(object):

    def __init__(self, letter, mass, name, intname, onetype):
        self.letter = letter
        self.ID     = []

        for i in range(len(nonglycine_expanded_sequence)):
            if self.letter == nonglycine_expanded_sequence[i]:
                self.ID.append(i)
        self.attrib = {'Mass':mass,'Name':str(name),'IntName':str(intname),"Type":str(onetype)}

#############################################
###               Set up XML              ###
#############################################

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

AAs = [] #List of Species objects representing each AA, in the same order as sites
for i in range (22):
    intname = 'Unbond ' + nonglycine_sites[i] + ' ' + nonglycine_sites[i]
    if ['NH', 'CH', 'CO'].count(nonglycine_sites[i]):
        intname = 'Backbone'
    AAs.append( Species( letter = nonglycine_sites[i], mass = PRIME20_masses.mass[nonglycine_sites[i]], name = nonglycine_names[i], intname = intname, onetype = 'Point') )

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
for species in AAs:
    #Add the species if it has any members
    if species.ID:
        temp = ET.SubElement( Genus, 'Species', attrib = species.attrib )
        temp = ET.SubElement( temp, 'IDRange', attrib = {'Type':'List'} )
        #Add its IDs
        [ ET.SubElement( temp, 'ID', attrib = { 'val' : str(i).replace("[","").replace("]","") } ) for i in species.ID ]

temp = ET.SubElement( Globals, 'Global', attrib = {'Type':'Cells','Name':'SchedulerNBList','NeighbourhoodRange':'7.400000000000e+00'})
ET.SubElement( temp, 'IDRange', attrib = {'Type':'All'})

#######################################################################
#               Create interaction types in XML format                #
#######################################################################

#iterate over all pairs of site types and write all their interactions
for n in range( len( nonglycine_sites ) ):
    site_type = nonglycine_sites[n]

    for m in range( n+1 ):
        site_type2 = nonglycine_sites[m]

        #Backbone-backbone interactions are handled by PRIME_BB interaction:
        if ['NH','CO','CH'].count(site_type2): 
            continue

        try:
            del attribset['Lambda']
            del attribset['WellDepth']
        except KeyError:
            pass

        pair             = site_type + ' ' + site_type2
        interaction_list = [] #At end of each inner loop, will be used to print which inters are defined for the pair.

        #################################
        #  Define unbonded interaction  #
        #################################

        attribset['Name']      = 'Unbond ' + pair
        attribset['Diameter']  = unbond_diam[pair]
        try:
            attribset['Lambda']    = unbond_lambda[pair]
            attribset['WellDepth'] = welldepth[pair]
            attribset['Type']      = 'SquareWell'
            unbondint[pair]        = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('SquareWell')
        except KeyError:
            attribset['Type'] = 'HardSphere'
            unbondint[pair]   = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('HardSphere')

        attribset = {'Elasticity':'1'}

        #######################################
        #  Define unbonded close interaction  #
        #######################################

        attribset['Name']     = 'Unbond close ' + pair
        attribset['Diameter'] = close_diam[pair]
        try:
            attribset['Lambda']    = close_lambda[pair]
            attribset['WellDepth'] = welldepth[pair]
            attribset['Type']      = 'SquareWell'
            closeunbondint[pair]   = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('SquareWell Close')
        except KeyError:
            attribset['Type']    = 'HardSphere'
            closeunbondint[pair] = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('HardSphere Close')

        attribset = {'Elasticity':'1'}

        ###############################
        #  Define bonded interaction  #
        ###############################

        try: #BB_SC only
            attribset['Name']     = 'Bond ' + pair
            attribset['Type']     = 'SquareBond'
            attribset['Diameter'] = BB_SC_bond_diam[pair]
            attribset['Lambda']   = BB_SC_bond_lambda[pair]
            bondint[pair]         = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('Bond')
        except KeyError as E:
            pass

        ###############################
        #  Define pseudo interaction  #
        ###############################

        try: #Only SC-NH and SC-CO use this.
            attribset['Name']     = 'Pseudo ' + pair
            attribset['Type']     = 'SquareBond'
            attribset['Diameter'] = BB_SC_pseudo_diam[pair]
            attribset['Lambda']   = BB_SC_pseudo_lambda[pair]
            pseudobondint[pair]   = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('Pseudo')
        except KeyError:
            pass

#######################################################################
#        Add each pair of species to its relevant interaction         #
#######################################################################

ET.SubElement( Interactions, 'Interaction', attrib = {'Type':'PRIME_BB', 'Name':'Backbone', 'Start':'0', 'End':str(len(expanded_sequence)-1) } )

#Iterate over every SC-SC and SC-BB pair of sites. ID2 will only be SCs.

for ID1, type1 in enumerate(expanded_sequence):
    for ID2 in sidechain_IDs:
        if ID2<=ID1 and ID1 in sidechain_IDs:
            continue

        type2=expanded_sequence[ID2]

        ThisInteraction = getInteraction(ID1, ID2, expanded_sequence, sequence)

        if debug:
            print "Assigning pair", ID1, expanded_sequence[ID1], ID2, expanded_sequence[ID2], "as", ThisInteraction.attrib['Name']

        try:
            ET.SubElement(ThisInteraction[0], 'IDPair', attrib = {'ID1':str(ID1), 'ID2':str(ID2)})
        except IndexError:
            ET.SubElement(ThisInteraction, 'IDPairRange', attrib = {'Type':'List'})
            ET.SubElement(ThisInteraction[0], 'IDPair', attrib = {'ID1':str(ID1), 'ID2':str(ID2)})

#Now remove interactions that don't apply to any pair.
for interaction_dict in [unbondint, closeunbondint, bondint, pseudobondint]:
    for key, interaction in interaction_dict.items():
        if len(interaction) != 0:
            Interactions.append( interaction )
        else:
            del interaction_dict[key]

######################
#  Create PSF files  #
######################

print "----------------------------------------------------"
print "WARNING PSF FILE-GENERATOR IS OUT OF DATE AND WRONG."
print "----------------------------------------------------"

psf_atoms_section = ""
psf_bonds_section = ""

#Backbone atoms
for i_res, res in enumerate(sequence):
    for i_local_atom, atom in enumerate(['NH', 'CH', 'CO']):
        i_atom = i_res*3 + i_local_atom
        psf_atoms_section += "{0: >8d} {1: <4} {2: <4d} {3: <4} {4: <4} {4: <4} {5: >10} {6: >13} {7: >11}\n".format(i_atom+1, str(0), i_res, res, atom, "0.000000", "0.0000", "0")

#SC atoms
for i_res, res in enumerate(nonglycine_expanded_sequence[n_bb_sites:]):
    i_atom = n_bb_sites + i_res
    psf_atoms_section += "{0: >8d} {1: <4} {2: <4d} {3: <4} {4: <4} {4: <4} {5: >10} {6: >13} {7: >11}\n".format(i_atom+1, str(0), i_res, res, res, "0.000000", "0.0000", "0")

#BB bonds
for i_bb_site in range(1,n_bb_sites):
    psf_bonds_section += "{0: >8d}{1: >8d}".format(i_bb_site, i_bb_site+1)

    if len(psf_bonds_section) - psf_bonds_section.rfind("\n") > 63:
        psf_bonds_section += "\n"

#SC bonds
for i_res, res in enumerate(sequence):
    if res != 'G':
        i_bb_site = i_res*3 + 2
        i_sc_site = n_bb_sites + 1 + i_res - sequence[:i_res].count('G')
        psf_bonds_section += "{0: >8d}{1: >8d}".format(i_bb_site, i_sc_site)

        if len(psf_bonds_section) - psf_bonds_section.rfind("\n") > 63:
            psf_bonds_section += "\n"

with open(psf_fn, 'w') as psf_file:
    psf_file.write("PSF\n\n\t1 !NTITLE\n REMARKS " + ''.join(sequence) + " STRUCTURE FILE\n REMARKS DATE: " + date + "\n\n")
    psf_file.write("{0: >8d}".format(n_bb_sites+n_sc_sites) + " !NATOM\n" + psf_atoms_section + "\n")
    psf_file.write("{0: >8d}".format(n_bb_sites-1+n_sc_sites) + " !NBOND\n" + psf_bonds_section + "\n\n")

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
    print subprocess.check_output(thermostat_command)
else:
    silent_stdout = subprocess.check_output(thermostat_command)

#Check config is valid with dynarun:
run_command = ['dynarun', '-c', '1000', '-o', xml_fn, xml_fn]
print "Running this command:", " ".join(run_command)
if debug:
    print subprocess.check_output(run_command)
else:
    silent_stdout = subprocess.check_output(run_command)

if debug:

    #Create trajectory file
    traj_command = ['dynamo2xyz', xml_fn]
    print "Running this command:", " ".join(traj_command)
    with open('traj.xyz', 'w') as trajfile:
        xyz = subprocess.check_output(traj_command)
        trajfile.write(xyz)

    convert_command = ["catdcd", "-o", dcd_temp_dir+"/dynamO_traj.dcd", "-xyz", "traj.xyz"]
    subprocess.check_output(convert_command)
    print "Running this command:", " ".join(convert_command)

    os.remove("traj.xyz")
