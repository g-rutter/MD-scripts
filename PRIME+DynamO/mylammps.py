from lammps import lammps
import numpy as np
import random

class mylammpsclass(lammps):

    def commands(self, iterable):
        """Run multiple commands, passed in an iterable."""

        for cmd in iterable:
            self.command(cmd)

        return 0

    def freshgroup(self, groupname, range):
        "Clear group and create again with new contents."

        self.command( "group " + groupname + " id 0" )
        self.command( "group " + groupname + " delete" )
        self.command( "group " + groupname + " " + range )

        return 0

    def minimize(self, etol=0.0, ftol=0.0, maxiter=1000, maxeval=10000):
        "Run a minimize."

        min_string = " ".join([ str(argument) for argument in (etol, ftol, maxiter, maxeval) ])
        self.command( "minimize "+min_string )

    def unfixes(self, iterable):
        """Remove list of fixes."""
        for fix in iterable:
            self.command( "unfix " + fix )

def create_groups(lmp, natoms, residue_atoms, sequence, i_res):
    """Create the standarded groups I need."""

    #Create group of last one defined residue
    last = natoms+len(residue_atoms)
    first = natoms+1
    lmp.freshgroup( "last1", "id "+str(first)+":"+str(last) )

    #Create group of last two defined residues
    last = natoms+len(residue_atoms)
    first = natoms-2
    if sequence[i_res-1] != 'G':
        first -= 1
    if first < 1:
        first = 1
    lmp.freshgroup( "last2", "id "+str(first)+":"+str(last) )

    #Create group of atoms involved in peptide bond
    NH = natoms+1
    CO = natoms
    if sequence[i_res-1] != 'G':
        CO -= 1
    lmp.freshgroup( "CO_NH", "id "+str(CO)+" "+str(NH) )

    #Create group for CH_NH pseudobond
    CH = natoms-1
    NH = natoms+1
    if sequence[i_res-1] != 'G':
        CH -= 1
    lmp.freshgroup( "CH_NH", "id "+str(CH)+" "+str(NH) )

    #Create group for CO_CH pseudobond
    CH = natoms+2
    CO = natoms
    if sequence[i_res-1] != 'G':
        CO -= 1
    lmp.freshgroup( "CO_CH", "id "+str(CH)+" "+str(CO) )

    #Create groups for addforce fix
    firstbbatom = 1
    lastbbatom  = natoms+3
    lmp.freshgroup( "firstbbatom", "id "+str(firstbbatom) )
    lmp.freshgroup( "lastbbatom",  "id "+str(lastbbatom)  )

def create_coords(coords, i_res, seq, res_space, prototype_positions):

    res = seq[i_res]
    if res == 'G':
        new_atoms = 3
    else:
        new_atoms = 4

    if i_res == 0:
        new_coords = prototype_positions[:new_atoms]

    else:
        prev_res = seq[i_res-1]
        if prev_res == 'G':
            last_new_atoms = 3
        else:
            last_new_atoms = 4

        i_last_NH = len(coords) - last_new_atoms

        #Place new NH at coords[i_last_NH]+res_space, place other atoms accordingly
        new_NH_pos = coords[i_last_NH] + [ res_space, 0.0, 0.0 ]
        shift = new_NH_pos - prototype_positions[0]

        new_size = coords.shape + np.array([new_atoms,0])
        new_coords = np.array(coords)
        new_coords.resize(new_size)
        new_coords[-new_atoms:] = prototype_positions[-new_atoms:] + shift

    return new_coords

