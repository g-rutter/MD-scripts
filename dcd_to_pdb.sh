#!/bin/bash

if [ $# -lt 3 ]; then
   echo "Usage: $0 [DCD_file] [PSF_file] [stride]"
   exit
fi

DCD_file=$1
DCD_filename=$(basename "$DCD_file")
DCD_filename_noext="${DCD_file%.*}"

PSF_file=$2
stride=$3

#Other settings
output_prefix=${DCD_filename_noext}"-"${stride}
PDB_traj=${output_prefix}".pdb"
temp=${output_prefix}".temp.pdb"

#Check DCD exists
if [ -f $DCD_file ]; then
   true
else
   echo "No such file $DCD_file. Exiting."
   exit
fi

##############################
#  Make PDB files from DCDs  #
##############################

echo "Converting DCD file to PDB file."
catdcd -o ${PDB_traj} -otype pdb -stride ${stride} -stype psf -s ${PSF_file} ${DCD_file}
echo "Using regex to fix the PDB trajectory file."
cat ${PDB_traj} | perl -pe 's{^[^A].*$}{$n=$n+1; "ENDMDL\nMODEL     $n"}e' | tail +2 | head --lines=-1 >$temp
cat $temp | perl -pe 's/(ATOM +\d+)    [41]/\1  N  /' | perl -pe 's/(ATOM +\d+)    2/\1  CA /' | perl -pe 's/(ATOM +\d+)    3/\1  C  /' | perl -pe 's/(ATOM +\d+) +\d+ +(\w)/\1  \2SC \2/' >${PDB_traj}
rm -f $temp
