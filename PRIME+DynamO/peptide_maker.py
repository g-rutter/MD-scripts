#!/usr/bin/python2.6

from lxml import etree as ET
import sys
import os
import time
import math
import random
import numpy as np

from mylammps import mylammpsclass
import mylammps

import PRIME20_unbonded
import PRIME20_bonded
import PRIME20_masses

def gauss():
    return str( random.gauss(0.0, 1.0) )

def getInteraction(ID1, ID2, expanded_sequence, n_bb_sites):
    """Checks if ID1 and ID2 are bonded or not, and assigns them an interaction from the correct dict of interactions."""
    pair1 = expanded_sequence[ID2] + ' ' + expanded_sequence[ID1]
    pair2 = expanded_sequence[ID1] + ' ' + expanded_sequence[ID2]

    #ID1 is always side-chain. ID1CH is the CH it's connected to.
    ID1CH = 3*( ID1 - n_bb_sites) + 1

    if ID2 < n_bb_sites:
        #ID2 is bb site

        if ID2 == ID1CH - 1: #SC pseudobonded to NH
            return pseudobondint[pair1]

        elif ID2 == ID1CH: #SC bonded to CH
            return bondint[pair1]

        elif ID2 == ID1CH + 1: #SC pseudobonded to CO
            return pseudobondint[pair1]

        #Unbonded if it gets to here. Check if it's a 'close' interaction.
        elif ( ID1CH - 3 ) < ID2 and ( ID1CH + 3 ) > ID2:
            try:
                return closeunbondint[pair2]
            except KeyError:
                return closeunbondint[pair1]

    #If both particles are sidechain, or if the above falls through, unbonded (and not 'close')
    try:
        return unbondint[pair2]
    except KeyError:
        return unbondint[pair1]

def joinstr(string, list_like):
    return string.join(str(item) for item in list_like)

#############################################
###       Define interaction dicts        ###
#############################################

#Define a dict of bonded and unbonded interactions, given as a pointer to the lxml SubElement variable.
unbondint      = {}
closeunbondint = {}
bondint        = {}
pseudobondint  = {}

attribset         = {'Elasticity':'1'}
unbond_diam       = dict(PRIME20_unbonded.SC_SC_diam)
unbond_welldepth  = dict(PRIME20_unbonded.SC_SC_well_depth)
unbond_welldiam   = dict(PRIME20_unbonded.SC_SC_well_diam)
bond_diam         = dict(PRIME20_bonded.BB_BB_diam)
pseudobond_diam   = dict(PRIME20_bonded.BB_pseudo)
pseudobond_lambda = '1.0450' #}http://www.sciencedirect.com/science/article/pii/S0022283611013672 Materials and Methods
bond_lambda       = '1.0450' #}
close_unbond_diam = dict(PRIME20_unbonded.close_diam)

unbond_diam.update(PRIME20_unbonded.BB_SC_diam)
unbond_diam.update(PRIME20_unbonded.BB_BB_diam)
bond_diam.update(PRIME20_bonded.BB_SC_diam)
unbond_welldepth.update(PRIME20_unbonded.BB_SC_well_depth)
unbond_welldiam.update(PRIME20_unbonded.BB_SC_well_diam)

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
    filename = 'Hall_peptide.xml'
try:
    temperature = sys.argv[2]
except:
    temperature = '1.0'

try:
    sys.argv[1] = [ letter.upper() for letter in sys.argv[1] ]
    if ( sys.argv[1] == list('N16N') ):
        sequence = list('AYHKKCGRYSYCWIPYDIERDRYDNGDKKC')
    else:
        sequence = list(sys.argv[1])
        assert [ sites.index(residue) for residue in sequence]
except:
    sys.exit('Run as ./peptide_maker.py (sequence) [temperature kT = 1.0] [filename = Hall_peptide].')

date                         = time.strftime('%X %x %Z')
box_size                     = len(sequence)*10.0
mid_box                      = box_size*0.5
box_pad                      = float(pseudobond_diam['CH CH'])
n_residues                    = len(sequence)
n_bb_sites                     = 3*n_residues
n_sc_sites                     = n_residues - sequence.count('G')
expanded_sequence            = ['NH','CH','CO']*n_residues + sequence #List of sites including G sites which dont really exist
nonglycine_expanded_sequence = filter(lambda a: a != 'G', expanded_sequence) #'Real' list, e.g. AGA gives NH,CH,CO,NH,CH,CO,NH,CH,CO,A,A
gap                          = 0.9

print 'Sequence:' , ''.join(sequence)
print 'File name:' , filename , '\n'

#############################################
###      Geometry of one amino acid       ###
#############################################

lmp = mylammpsclass()

lmp.command( "units real" )
lmp.command( "dimension 3" )
lmp.command( "boundary m m m" )
lmp.command( "atom_style molecular" )
lmp.command( "atom_modify sort 0 0" )

lmp.command( "read_data box.lmp" )
lmp.file("LAMMPS_PRIME20.params")

lmp.command( "pair_style lj/cut 10.0" )
lmp.command( "pair_coeff * * 0.0 10.0" )
lmp.command( "special_bonds lj 0 1 1" )

lmp.command( "fix NVT all nvt temp 300 0.1 100.0" )

atom_tally = 0

for i_res, res in enumerate(sequence):
    natoms = lmp.get_natoms()

    lmp_coords = np.array(lmp.gather_atoms( "x", 1, 3))
    lmp_coords = np.reshape(lmp_coords, [natoms, 3])
    new_coords = mylammps.create_coords(lmp_coords, i_res, res)

    #Add next residue
    residue_atoms = [ 'NH', 'CH', 'CO']
    if res != 'G':
        residue_atoms += res

    #Check atoms have been added correctly
    if atom_tally != natoms:
        print "ERROR: Internal atom tally doesn't match LAMMPS' natoms."
        print "Atom tally:", atom_tally
        print "LAMMPS natoms:", natoms
        exit()

    atom_tally += len(residue_atoms)

    #Create atoms for next residue
    for i_atom, atom in enumerate(residue_atoms):
        coordstr = joinstr(" ", new_coords[natoms+i_atom])
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

    pullstrength = 100.0*i_res
    forcecmds.append( "fix  negativepull  firstbbatom addforce -{0:f} 0.0 0.0". format(pullstrength) )
    forcecmds.append( "fix  positivepull  lastbbatom  addforce  {0:f} 0.0 0.0". format(pullstrength) )
    force_fix_names = ["negativepull", "positivepull"]

    lmp.commands(bondcmds+forcecmds)

    lmp.command( "dump dcd all dcd 20 " + str(i_res) + "test.dcd" )

    #Run
    lmp.command( "run 1" )
    lmp.minimize()

    lmp.command( "run 2000" ) 
    lmp.command( "undump dcd" )
    lmp.minimize()

    lmp.unfixes(bond_fix_names+force_fix_names)

natoms = lmp.get_natoms()

lmp_coords = np.array(lmp.gather_atoms( "x", 1, 3))
lmp_coords = np.reshape(lmp_coords, [natoms, 3])

lmp.close()

#############################################
###          Obtain all positions         ###
#############################################

#This prog mainly puts all the SC sites at the end, while the LAMMPS portion
#includes them at the end of each residue, so some translation is required.

i_sc = i_lammps_atom = 0

bb_positions = np.zeros([n_bb_sites, 3])
sc_positions = np.zeros([n_residues, 3]) #sparsely populated

for i_res, res in enumerate( sequence ):
    #Add next residue
    residue_atoms = [ 'NH', 'CH', 'CO']

    for i_atom in [0,1,2]:
        bb_positions[i_res*3 + i_atom] = lmp_coords[i_lammps_atom]
        i_lammps_atom += 1

    if res != 'G':
        sc_positions[i_res] = lmp_coords[i_lammps_atom]
        i_lammps_atom += 1

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

[ ET.SubElement( ParticleData, 'Pt', attrib = {'ID' : str(ID) } ) for ID in range( n_bb_sites + n_sc_sites ) ]
#Give positions and momenta to the backbone sites:
for i in range( n_bb_sites ):
    ET.SubElement( ParticleData[i], 'P', attrib = {"x":str( bb_positions[i][0] ), "y":str( bb_positions[i][1] ), "z":str( bb_positions[i][2] )} )
    ET.SubElement( ParticleData[i], 'V', attrib = {"x":gauss(), "y":gauss(), "z":gauss()} )

#Give positions and momenta to the side chain sites:
j=0
for i in range( len(sequence) ):
    if sequence[i] != 'G':
        ET.SubElement( ParticleData[j + n_bb_sites], 'P', attrib = {"x":str( sc_positions[i][0] ), "y":str( sc_positions[i][1] ), "z":str( sc_positions[i][2] )} )
        ET.SubElement( ParticleData[j + n_bb_sites], 'V', attrib = {"x":gauss(), "y":gauss(), "z":gauss()} )
        j += 1

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

temp = ET.SubElement( Globals, 'Global', attrib = {'Type':'Cells','Name':'SchedulerNBList','NeighbourhoodRange':'5.000000000000e+00'})
ET.SubElement( temp, 'IDRange', attrib = {'Type':'All'})

#######################################################################
#               Create interaction types in XML format                #
#######################################################################

for n in range( len( nonglycine_sites ) ):
    site_type = nonglycine_sites[n]

    try: #Don't bother writing interactions for site types that don't exist.
        expanded_sequence.index(site_type)
    except ValueError:
        continue

    for m in range( n+1 ):
        site_type2 = nonglycine_sites[m]

        if ['NH','CO','CH'].count(site_type2): #Backbone-backbone interactions are handled by PRIME_BB interaction.
            continue

        try: #Don't bother writing interactions for site types that don't exist.
            expanded_sequence.index(site_type2)
        except ValueError:
            continue

        try:
            del attribset['Lambda']
            del attribset['WellDepth']
        except KeyError:
            pass

        pair             = site_type + ' ' + site_type2
        interaction_list = [] #At end of each inner loop, will be used to print which inters are defined for the pair.

        attribset['Name']      = 'Unbond ' + pair
        attribset['Diameter']  = unbond_diam[pair]
        try:
            attribset['Lambda']    = str ( float( unbond_welldiam[pair] ) / float( unbond_diam[pair] ) )
            attribset['WellDepth'] = unbond_welldepth[pair]
            attribset['Type']      = 'SquareWell'
            unbondint[pair]        = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('SquareWell')
        except KeyError:
            attribset['Type'] = 'HardSphere'
            unbondint[pair]   = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('HardSphere')

        try:
            del attribset['Lambda']
            del attribset['WellDepth']
        except KeyError:
            pass

        attribset['Name']     = 'Unbond close ' + pair
        attribset['Diameter'] = close_unbond_diam[pair]
        try:
            attribset['Lambda']    = str ( float( unbond_welldiam[pair] )/float( close_unbond_diam[pair] ) )
            attribset['WellDepth'] = unbond_welldepth[pair]
            attribset['Type']      = 'SquareWell'
            closeunbondint[pair]   = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('SquareWell Close')
        except KeyError:
            attribset['Type']    = 'HardSphere'
            closeunbondint[pair] = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('HardSphere Close')

        try:
            del attribset['Lambda']
            del attribset['WellDepth']
        except KeyError:
            pass

        try: #No bonded interaction exists for many SC cases.
            attribset['Name']     = 'Bond ' + pair
            attribset['Type']     = 'SquareBond'
            attribset['Diameter'] = bond_diam[pair]
            attribset['Lambda']   = bond_lambda
            bondint[pair]         = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('Bond')
        except KeyError:
            pass

        try: #No pseudo interaction defined for many pairs.
            attribset['Name']     = 'Pseudo ' + pair
            attribset['Type']     = 'SquareBond'
            attribset['Diameter'] = pseudobond_diam[pair]
            attribset['Lambda']   = pseudobond_lambda
            pseudobondint[pair]   = ET.Element ( 'Interaction', attrib = attribset )
            interaction_list.append('Pseudo')
        except KeyError:
            pass

#######################################################################
#        Add each pair of species to its relevant interaction         #
#######################################################################

ET.SubElement( Interactions, 'Interaction', attrib = {'Type':'PRIME_BB', 'Name':'Backbone', 'Start':'0', 'End':str(n_bb_sites-1) } )

#Iterate over every SC-SC and SC-BB pair of sites. ID1 will only be SCs.
#ID1 and ID2 include glycine as if its a real site, so IDx_flat corrects this.

for ID1 in range( n_bb_sites, len(expanded_sequence) ):

    if expanded_sequence[ID1] == 'G':
        continue

    ID1_flat = ID1 - expanded_sequence[:ID1].count('G')

    for ID2 in range( ID1 + 1):

        if expanded_sequence[ID2] == 'G':
            continue

        ID2_flat= ID2 - expanded_sequence[:ID2].count('G')

        #getInteraction assumes no glycine, i.e. the last SC atom will have same id in GA as in AA.
        ThisInteraction = getInteraction(ID1, ID2, expanded_sequence, n_bb_sites)
        #For debug:
        #print ("Assigning pair", ID1, expanded_sequence[ID1], ID2, expanded_sequence[ID2], "as", ThisInteraction.attrib['Name'])

        try:
            ET.SubElement(ThisInteraction[0], 'IDPair', attrib = {'ID1':str(ID1_flat), 'ID2':str(ID2_flat)})
        except IndexError:
            ET.SubElement(ThisInteraction, 'IDPairRange', attrib = {'Type':'List'})
            ET.SubElement(ThisInteraction[0], 'IDPair', attrib = {'ID1':str(ID1_flat), 'ID2':str(ID2_flat)})

        #print ID1_flat, ID2_flat

#Now remove interactions that don't apply to any pair.
for interaction_dict in [unbondint, closeunbondint, bondint, pseudobondint]:
    for key, interaction in interaction_dict.items():
        if len(interaction) != 0:
            Interactions.append( interaction )
        else:
            del interaction_dict[key]

#############################################
###              Write file               ###
#############################################

input_file = open(filename, 'w')
input_file.write('<!-- DynamO input file contains the PRIME20 model of the sequence: ' + ''.join(sequence) + '. -->\n')
input_file.write('<!-- Created on ' +date + '. -->\n')
[ input_file.write(ET.tostring(DynamOconfig, pretty_print=True)) ]
input_file.close()

#Add thermostat and rescale via dynamod:
thermostat_command = 'dynamod -T ' + temperature + ' -r ' + temperature + ' -o ' + filename + ' -Z ' + filename + " --round"
print "Running this command:", thermostat_command
os.system(thermostat_command)
