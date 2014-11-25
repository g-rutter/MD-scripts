#!/bin/sh

##############
#  Settings  #
##############

xyz=temp.xyz
pdb=out.pdb
temp=temp.temp
max_args=1000

##################
#  Get psf file  #
##################

if [ -z "$1" ]; then
    psf=(*.psf)
else
    psf=$1
fi

echo "PSF = $psf."
echo "If wrong pass it as the first argument."

###################################
#  Collect and order input files  #
###################################

all_files=(Snapshot.*.xml.bz2)
snapshot_files=${all_files[@]//*output*}

ordered_snapshot_files=(`for file in ${snapshot_files[@]}; do echo "$file"; done\
    | sort -n --field-separator=. -k2`)

n_snapshot_files=${#ordered_snapshot_files[*]}

echo $n_snapshot_files files found.

###################
#  Make xyz file  #
###################

echo "Making XYZ file."

start_nos=`seq 0 $max_args $n_snapshot_files`

rm -f $xyz

for start_no in ${start_nos[@]}; do
    dynamo2xyz ${ordered_snapshot_files[@]:$start_no:$max_args} >>$xyz
done

echo ""

###################
#  Make pdb file  #
###################

echo "Making PDB file."
catdcd -o ${pdb} -otype pdb -stype psf -s $psf -xyz $xyz
echo ""

echo "Using regex to fix the PDB trajectory file."
cat ${pdb} | perl -pe 's{^[^A].*$}{$n=$n+1; "ENDMDL\nMODEL     $n"}e' | sed 1d | sed '$d' >$temp
cat $temp | perl -pe 's/(ATOM +\d+)    [41]/\1  N  /' | perl -pe 's/(ATOM +\d+)    2/\1  CA /' | perl -pe 's/(ATOM +\d+)    3/\1  C  /' | perl -pe 's/(ATOM +\d+) +\d+ +(\w)/\1  \2SC \2/' >${pdb}
rm -f $temp $xyz
echo ""

echo "Done. Check: $pdb"
