#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 [PDB file] [N_chains]"
    exit
fi

################################
#  Settings from command-line  #
################################

PDB_file=$1
N=$2

###################
#  Auto-settings  #
###################


PDB_filename_noext="${PDB_file%.*}"
output_prefix=${PDB_filename_noext}"-splitchains"
frames_in=`catdcd -num -pdb $PDB_file | ack "^Total" | cut -d " " -f 3`

echo ${frames_in} frames

###########################
#  Get M atoms in system  #
###########################

ENDMDL_line_nr=`awk '/ENDMDL/{print NR; exit}' $PDB_file`
MDL1_last_line=`head --lines=$((ENDMDL_line_nr-1)) $PDB_file | tail -1`
M_atoms=`echo $MDL1_last_line | cut -d ' ' -f 2`
atoms_per_chain=$((M_atoms/N))

echo -e "Detected:\n$M_atoms atoms total, making up\n$N chains with\n$atoms_per_chain atoms each.\n"

###########
#  Split  #
###########

echo "Writing new file..."

awk '

BEGIN{
    model_i  = 2
    mol_last = 0
    atoms_per_chain = '$atoms_per_chain'

    printf "MODEL     1\n"
}

/ATOM/{
    mol_float=($2-1)/atoms_per_chain
    mol_int=int(mol_float)

    if (mol_int != mol_last){
        printf "ENDMDL\n"
        printf "MODEL     %-i\n", model_i
        model_i ++
    }

    atom_id = $2 - (atoms_per_chain*mol_int)
    model_id = 0

    printf("ATOM  %5i  %-3s %-3s %-2i %-6i%8.3f%8.3f%8.3f%6.2f%6.2f      %-4i  \n", atom_id, $3, $4, model_id, $6, $7, $8, $9, $10, $11, $12);

    mol_last = mol_int
}

END{
    printf "ENDMDL\n"
}


' $PDB_file >$output_prefix.pdb

echo "Done."
