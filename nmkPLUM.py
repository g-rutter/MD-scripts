#!/usr/bin/python

#Back-bone atom type-number | 1        | 2            | 3       | 4
#Corresponding atom         | Nitrogen | Carbon alpha | Carbon' | Proline's nitrogen

#Side-chain type-number | 5            | 6             | 7                | 8                | 9                | 10
#Residue represented    | Alanine(A)   | Cysteine(C)   | Aspartic acid(D) | Glutamic acid(E) | Phenylalanine(F) | Histidine(H)

#Side-chain type-number | 11            | 12        | 13           | 14            | 15             | 16          | 17
#Residue represented    | Isoleucine(I) | Lysine(K) | Leucine(L)   | Methionine(M) | Asparagine(N)  | Proline(P)  | Glutamine(Q)

#Side-chain type-number | 18            | 19        | 20           | 21            | 22             | 23
#Residue represented    | Arginine(R)   | Serine(S) | Threonine(T) | Valine(V)     | Tryptophan (W) | Tyrosine(Y)

import numpy, random, sys, time, math
import re

##############
#  Settings  #
##############
# Size units are Angstroms = 0.1 nanometers

#This is the padding around the smallest box containing all atoms in the initial configuration.
box_pad  = 15.0

#Displacement in x, y and z between the first N site of each residue and the next.
peptide_pad = 15.0

##################################
#  Read in command-line options  #
##################################

SC         = list('ACDEFGHIKLMNPQRSTVWY')
nongly_SC  = list('ACDEFHIKLMNPQRSTVWY')

outputs = { 'stdout' : sys.stdout, 'stderr' : sys.stderr}
help = "Run as: " + sys.argv[0] + " [sequence] [number] {filename}"

try:
   N_proteins=int(sys.argv[2])
except IndexError, ValueError:
   exit(help)

#set and check sequence
try:
    sys.argv[1] = sys.argv[1].upper()
    if ( sys.argv[1] == 'S1' ):
        sequence = list('PPPWLPYMPPWS')
    elif ( sys.argv[1] == 'N16N' ):
        sequence = list('AYHKKCGRYSYCWIPYDIERDRYDNGDKKC')
    elif ( sys.argv[1] == 'N16NN' ):
        sequence = list('AYHKKCGRYSYCWIPYNIQRNRYNNGNKKC')
    elif ( sys.argv[1][:5] == 'ALPHA' ):
        sequence = list('ACDEFGHIKLMNPQRSTVWY')
    elif (sys.argv[1] == '2A3D'):
        sequence = list('MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN')
    else:
        sequence = list(sys.argv[1])
        valid_input = SC + [ str(number) for number in range(0,10) ]
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
    exit(help)

try:
   lmp_file = outputs[ sys.argv[3] ]
   psf_file = outputs[ sys.argv[3] ]
except KeyError:
   filename = sys.argv[3]
   lmp_file = open(filename+".lmp", 'w')
   psf_file = open(filename+".psf", 'w')
except IndexError:
   lmp_file = open(sys.argv[1] + "-" + sys.argv[2] + ".lmp", 'w')
   psf_file = open(sys.argv[1] + "-" + sys.argv[2] + ".psf", 'w')

######################
#  Define functions  #
######################

def randomVector(length):

    rv=numpy.array([random.uniform(-1.0,1.0),random.uniform(-1.0,1.0),random.uniform(-1.0,1.0)])
    norm=length/numpy.sqrt(numpy.dot(rv,rv))
    rv*=norm
    return rv

def makeCoords(sequence,first_position):

    #Create array
    positions=numpy.zeros((N_protein_sites,3))

    #Place N
    positions[0,0]=first_position[0]
    positions[0,1]=first_position[1]
    positions[0,2]=first_position[2]

    #Place other atoms
    i0=0
    iLast=0
    iRes=0

    dr=randomVector(1.455)

    positions[i0+1,:]=positions[i0,:]+dr

    dr=randomVector(1.51)
    positions[i0+2,:]=positions[i0+1,:]+dr

    ## side chain (for non-glycine residues)
    try:
        at=str(nongly_SC.index(sequence[iRes])+5)
        dr=randomVector(1.53)

        positions[i0+3,:]=positions[i0+1,:]+dr

    except:
        pass

    
    if sequence[iRes]=='G':
        i0+=3
    else:
        i0+=4

    for iRes in range(1,len(sequence)):
        ## place N
        if sequence[iRes-1]=='G':
            iLast=i0-1
        else:
            iLast=i0-2

        dr=randomVector(1.325)

        positions[i0,:]=positions[iLast,:]+dr

        ##  place Ca
        dr=randomVector(1.455)
        positions[i0+1,:]=positions[i0,:]+dr


        ## place C'
        dr=randomVector(1.51)
        positions[i0+2,:]=positions[i0+1,:]+dr

        ## side chain (for non-glycine residues)
        try:
            dr=randomVector(1.53)
            positions[i0+3,:]=positions[i0+1,:]+dr

        except:
            pass

        if sequence[iRes]=='G':
            i0+=3
        else:
            i0+=4

    return positions

N_protein_sites=0
bond_count_per=0
angle_count_per=0
dihedral_count=0

## initialise a few bits and bobs

date = time.strftime('%X %x %Z')

## get number of beads in protein
N_protein_sites=0
for residue in sequence:
    if residue=='G':
        N_protein_sites+=3
    else:
        N_protein_sites+=4

#####################
#  Place all sites  #
#####################

points_per_dimension = int( math.ceil(N_proteins**(1.0/3)) )
atoms=numpy.zeros((N_proteins, N_protein_sites, 3))
first_position = [0.0,0.0,0.0]

box_core = peptide_pad * points_per_dimension
box_bounds = [ -box_core/2, box_core/2 ]

for iPro in range(N_proteins):

    xpoint = iPro % points_per_dimension
    ypoint = iPro/(points_per_dimension**2)
    zpoint = (iPro/points_per_dimension)%points_per_dimension

    try:
       first_position[0] = box_bounds[0] + ( box_core*(float(xpoint)/(points_per_dimension-1)) )
       first_position[1] = box_bounds[0] + ( box_core*(float(ypoint)/(points_per_dimension-1)) )
       first_position[2] = box_bounds[0] + ( box_core*(float(zpoint)/(points_per_dimension-1)) )
    except ZeroDivisionError:
       first_position[0] = box_bounds[0]
       first_position[1] = box_bounds[0]
       first_position[2] = box_bounds[0]

    atoms[iPro]=makeCoords(sequence,first_position)

#expand inner_bounds to include every site
box_bounds[0] = -max(numpy.amax(atoms), -numpy.amin(atoms))
box_bounds[1] = -box_bounds[0]

atomType=numpy.zeros((N_protein_sites))
atomResidue=[]
atomResNumber=[]
bonds=[]
bondType=[]
angles=[]
angleType=[]
dihedrals=[]
dihedralType=[]
i0=0
iLast=0
iRes=0

## first residue

## place N
if sequence[iRes]=='P':
    atomType[i0]=4
else:
    atomType[i0]=1
atomResidue.append(sequence[iRes])
atomResNumber.append(1)

atomType[i0+1]=2
atomResidue.append(sequence[iRes])
atomResNumber.append(1)

atomType[i0+2]=3
atomResidue.append(sequence[iRes])
atomResNumber.append(1)

## side chain (for non-glycine residues)
try:
    at=str(nongly_SC.index(sequence[iRes])+5)

    atomType[i0+3]=at
    atomResidue.append(sequence[iRes])
    atomResNumber.append(1)
except:
    pass

## create bonds
bonds.append((i0+1,i0+2))
bondType.append(1)
bonds.append((i0+2,i0+3))
bondType.append(2)
if sequence[iRes] in nongly_SC:
    bonds.append((i0+2,i0+4))
    bondType.append(4)

## create angles
angles.append((i0+1,i0+2,i0+3))
angleType.append(1)
if sequence[iRes] in nongly_SC:
    angles.append((i0+1,i0+2,i0+4))
    angleType.append(4)
    angles.append((i0+3,i0+2,i0+4))
    angleType.append(5)

if sequence[iRes] in nongly_SC:
    dihedrals.append((i0+1,i0+2,i0+3,i0+4))
    dihedralType.append(5)

    
if sequence[iRes]=='G':
    i0+=3
else:
    i0+=4

for iRes in range(1,len(sequence)):
    ## place N
    if sequence[iRes-1]=='G':
        iLast=i0-1
    else:
        iLast=i0-2

    if sequence[iRes]=='P':
        atomType[i0]=4
    else:
        atomType[i0]=1
    atomResidue.append(sequence[iRes])
    atomResNumber.append(iRes+1)

    ##  place Ca
    atomType[i0+1]=2
    atomResidue.append(sequence[iRes])
    atomResNumber.append(iRes+1)

    ## place C'
    atomType[i0+2]=3
    atomResidue.append(sequence[iRes])
    atomResNumber.append(iRes+1)

    ## side chain (for non-glycine residues)
    try:
        at=str(nongly_SC.index(sequence[iRes])+5)

        atomType[i0+3]=at
        atomResidue.append(sequence[iRes])
        atomResNumber.append(iRes+1)

    except:
        pass

    bonds.append((iLast+1,i0+1))
    bondType.append(3)
    bonds.append((i0+1,i0+2))
    bondType.append(1)
    bonds.append((i0+2,i0+3))
    bondType.append(2)
    if sequence[iRes] in nongly_SC:
        bonds.append((i0+2,i0+4))
        bondType.append(4)

    ## create angles
    angles.append((iLast,iLast+1,i0+1))
    angleType.append(2)
    angles.append((iLast+1,i0+1,i0+2))
    angleType.append(3)
    
    angles.append((i0+1,i0+2,i0+3))
    angleType.append(1)
    if sequence[iRes] in nongly_SC:
        angles.append((i0+1,i0+2,i0+4))
        angleType.append(4)
        angles.append((i0+3,i0+2,i0+4))
        angleType.append(5)

    dihedrals.append((iLast-1,iLast,iLast+1,i0+1))
    dihedralType.append(1)
    dihedrals.append((iLast,iLast+1,i0+1,i0+2))
    if sequence[iRes]=='P':
        dihedralType.append(3)
    else:
        dihedralType.append(2)
    dihedrals.append((iLast+1,i0+1,i0+2,i0+3))
    dihedralType.append(4)
    if sequence[iRes] in nongly_SC:
        dihedrals.append((i0+1,i0+2,i0+3,i0+4))
        dihedralType.append(5)

    if sequence[iRes]=='G':
        i0+=3
    else:
        i0+=4

bond_count_per=len(bonds)
angle_count_per=len(angles)
dihedral_count=0
improper_count=0

###############
#  Write psf  #
###############

## write out psf file
psf_atoms_section=""
psf_bonds_section  = ""
psf_dihedrals_section = ""
psf_impropers_section = ""
lmp_dihedrals_section = ""

nRes=len(sequence)
for iPro in range(N_proteins):

    atoms_tot = iPro*N_protein_sites
    for i in range(N_protein_sites):
        psf_atoms_section += "{0: >8d} {1: <4} {2: <4d} {3: <4} {4: <4} {4: <4} {5: >10} {6: >13} {7: >11}\n".format(atoms_tot+i+1, str(iPro), iPro*nRes+atomResNumber[i], atomResidue[i], int(atomType[i]), "0.000000", "0.0000", "0")

    for bond in bonds:
        psf_bonds_section += "{0: >8d}{1: >8d}".format(atoms_tot+bond[0],atoms_tot+bond[1])
        if len(psf_bonds_section) - psf_bonds_section.rfind("\n") > 63:
            psf_bonds_section += "\n"

    #Adding dihedrals to psf file
    for i in range ( len ( dihedrals ) ):

        if dihedralType[i] is 5:
            improper_count += 1
            psf_impropers_section += "{0: >8d}{1: >8d}{2: >8d}{3: >8d}".format(atoms_tot+dihedrals[i][0],atoms_tot+dihedrals[i][1],atoms_tot+dihedrals[i][2],atoms_tot+dihedrals[i][3])
            if len(psf_impropers_section) - psf_impropers_section.rfind("\n") > 63:
                  psf_impropers_section += "\n"
        else:
            dihedral_count += 1
            psf_dihedrals_section += "{0: >8d}{1: >8d}{2: >8d}{3: >8d}".format(atoms_tot+dihedrals[i][0],atoms_tot+dihedrals[i][1],atoms_tot+dihedrals[i][2],atoms_tot+dihedrals[i][3])
            if len(psf_dihedrals_section) - psf_dihedrals_section.rfind("\n") > 63:
                  psf_dihedrals_section += "\n"

        lmp_dihedrals_section += "%6d %6d %6d %6d %6d %6d\n" % (dihedral_count+improper_count,dihedralType[i],atoms_tot+dihedrals[i][0],atoms_tot+dihedrals[i][1],atoms_tot+dihedrals[i][2],atoms_tot+dihedrals[i][3])

psf_file.write("PSF\n\n\t1 !NTITLE\n REMARKS " + ''.join(sequence) + " STRUCTURE FILE\n REMARKS DATE: " + date + "\n\n")
psf_file.write("{0: >8d}".format(N_proteins*N_protein_sites) + " !NATOM\n" + psf_atoms_section + "\n")
psf_file.write("{0: >8d}".format(N_proteins*bond_count_per) + " !NBOND\n" + psf_bonds_section + "\n\n")
psf_file.write("{0: >8d}".format(dihedral_count) + " !NPHI\n" + psf_dihedrals_section + "\n\n")
psf_file.write("{0: >8d}".format(improper_count) + " !NIMPHI\n" + psf_impropers_section)

psf_file.close()

###############
#  Write lmp  #
###############

print >> lmp_file, "PLUM CG protein" 
print >> lmp_file
print >> lmp_file, "%6d atoms" % (N_proteins*N_protein_sites)
print >> lmp_file, "%6d bonds" % (N_proteins*bond_count_per) 
print >> lmp_file, "%6d angles" % (N_proteins*angle_count_per)
print >> lmp_file, "%6d dihedrals" % (dihedral_count+improper_count)
print >> lmp_file
print >> lmp_file, "%6d atom types" % (23)
print >> lmp_file, "%6d bond types" % (4)
print >> lmp_file, "%6d angle types" % (5)
print >> lmp_file, "%6d dihedral types" % (5)
print >> lmp_file
print >> lmp_file, "%8.3f %8.3f xlo xhi" % (box_bounds[0]-box_pad, box_bounds[1]+box_pad)
print >> lmp_file, "%8.3f %8.3f ylo yhi" % (box_bounds[0]-box_pad, box_bounds[1]+box_pad)
print >> lmp_file, "%8.3f %8.3f zlo zhi" % (box_bounds[0]-box_pad, box_bounds[1]+box_pad)
print >> lmp_file
print >> lmp_file, "Atoms"
print >> lmp_file

for iPro in range(N_proteins):
    for i in range(N_protein_sites):
        xx=atoms[iPro][i,0]
        yy=atoms[iPro][i,1]
        zz=atoms[iPro][i,2]
        print >> lmp_file, "%6d %6d %6d %8.3f %8.3f %8.3f" % (iPro*N_protein_sites+i+1,iPro+1,atomType[i],xx,yy,zz)

print >> lmp_file
print >> lmp_file, "Bonds"
print >> lmp_file
for iPro in range(N_proteins):
    for i,(bt,b) in enumerate(zip(bondType,bonds)):

        print >> lmp_file, "%6d "*4 % (iPro*bond_count_per+i+1,bt,iPro*N_protein_sites+b[0],iPro*N_protein_sites+b[1])

print >> lmp_file
print >> lmp_file, "Angles"
print >> lmp_file
for iPro in range(N_proteins):
    for i,(at,a) in enumerate(zip(angleType,angles)):
        print >> lmp_file, "%6d "*5 % (iPro*angle_count_per+i+1,at,iPro*N_protein_sites+a[0],iPro*N_protein_sites+a[1],iPro*N_protein_sites+a[2])

print >> lmp_file
print >> lmp_file, "Dihedrals"
print >> lmp_file, lmp_dihedrals_section
print >> lmp_file
