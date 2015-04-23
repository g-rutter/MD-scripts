#!/bin/sh

if ! [ $# -eq 3 ]; then
    echo "Usage: $0 [target_cluster] [clust_id_file] [PDB]"
    exit
fi

target_cluster=$1
clust_id_file=$2
PDB=$3

tempdir="temp"
temp_pdb="temp.pdb"
out_pdb="C$target_cluster.pdb"

frames=(`awk '/ +[0-9]+ +'$target_cluster' *$/{print $1}' $clust_id_file`)

mkdir -p $tempdir

i=0

echo Grabbing ${#frames[*]} frames.

for frame in ${frames[*]}; do
    catdcd -stype pdb -s $PDB -otype pdb -o $tempdir/C$target_cluster-$i.pdb -first $frame -last $frame -pdb $PDB >/dev/null
    i=$((i+1))
done

echo Stitching together.

catdcd -stype pdb -s $PDB -otype pdb -o C$target_cluster.pdb -pdb $tempdir/C$target_cluster-*.pdb >/dev/null

rm $tempdir/C$target_cluster-*.pdb
rmdir --ignore-fail-on-non-empty $tempdir

echo "Using regex to fix up PDB."
cat ${out_pdb} | perl -pe 's{^[^A].*$}{$n=$n+1; "ENDMDL\nMODEL     $n"}e' | sed 1d | sed '$d' >$temp_pdb
cat $temp_pdb | perl -pe 's/(ATOM +\d+)    [41]/\1  N  /' | perl -pe 's/(ATOM +\d+)    2/\1  CA /' | perl -pe 's/(ATOM +\d+)    3/\1  C  /' | perl -pe 's/(ATOM +\d+) +\d+ +(\w)/\1  \2SC \2/' >${out_pdb}
rm -f $temp_pdb
