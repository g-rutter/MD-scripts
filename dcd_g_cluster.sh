#!/bin/sh

#################
#  Parse input  #
#################

if [ $# -lt 3 ]; then
   echo "Usage: $0 [DCD_file] [PSF_file] [stride] {clust_cut=0.4}"
   exit
elif [ $# -lt 4 ]; then
   clust_cut=0.4
else
   clust_cut=$4
fi

DCD_file=$1
DCD_filename=$(basename "$DCD_file")
DCD_filename_noext="${DCD_file%.*}"

PSF_file=$2
stride=$3

#Other settings
IDX_file="index.ndx"
output_prefix=${DCD_filename_noext}"-"${stride}
PDB_traj=${output_prefix}".pdb"
temp=${output_prefix}".temp.pdb"

#Check DCD exists
if [ -f $DCD_file ]; then
   continue
else
   echo "No such file $DCD_file. Exiting."
   exit
fi

#Check index file exists
if ! [ -f $IDX_file ]; then
   if ! [ -f ../$IDX_file ]; then
      echo "No index file named $IDX_file found. Exiting."
      exit
   fi
   IDX_file=../$IDX_file
fi

echo Settings:
echo DCD file: $DCD_file
echo PSF file: $PSF_file
echo Stride for reading from DCD file: $stride
echo Prefix for output files: $output_prefix
echo Cutoff RMSD for g_cluster: $clust_cut
echo Index file for g_cluster: $IDX_file
sleep 2s

##############################
#  Make PDB files from DCDs  #
##############################

if ! [ -f ${PDB_traj} ]; then
   echo "Converting DCD file to PDB file."
   catdcd -o ${PDB_traj} -otype pdb -stride ${stride} -stype psf -s ${PSF_file} ${DCD_file}
   echo "Using regex to fix the PDB trajectory file."
   cat ${PDB_traj} | perl -pe 's{^[^A].*$}{$n=$n+1; "ENDMDL\nMODEL     $n"}e' | tail +2 | head --lines=-1 >$temp
   cat $temp | perl -pe 's/(ATOM +\d+)    [41]/\1  N  /' | perl -pe 's/(ATOM +\d+)    2/\1  CA /' | perl -pe 's/(ATOM +\d+)    3/\1  C  /' | perl -pe 's/(ATOM +\d+) +\d+ +(\w)/\1  \2SC \2/' >${PDB_traj}
   rm -f $temp
else
   echo "${PDB_traj} already exists. Please delete to recalculate."
fi

###############
#  g_cluster  #
###############

if ! [ -f ${output_prefix}_clusters.pdb -a -f ${output_prefix}_clust_id.xvg ]; then
   echo "Clustering."
   expect -c "spawn g_cluster -s ${PDB_traj} -f ${PDB_traj} -method gromos -n ${IDX_file} -cl ${output_prefix}_clusters.pdb -cutoff ${clust_cut} -clid ${output_prefix}_clust_id.xvg -o ${output_prefix}_rmsd-clust.xpm -g ${output_prefix}_clust.log -dist ${output_prefix}_rmsd-dist.xvg

              expect \"Select group for least squares fit and RMSD calculation:\"
              expect \"3*Backbone\"
              expect \"Select a group: $\"
              send \"3\r\"

              expect \"Select group for output:\"
              expect \"4*All\"
              expect \"Select a group: $\"
              send \"4\r\"

              interact
             "
else
   echo "Clustering results already exist. Please delete ${output_prefix}_clusters.pdb to recalculate."
fi

###########################
#  Cluster count vs time  #
###########################

CLIDs=(`tail +16 ${output_prefix}_clust_id.xvg | sed 's/^ *[0-9]\+ \+//'`)
Unique_CLIDs=()
N_uniq_CLIDs=0
CLID_count=()

if ! [ -f ${output_prefix}_nclusters.xvg ]; then

   echo "Calculating number of clusters over time in file ${output_prefix}_nclusters.xvg"
   echo "Calculating population of clusters in file ${output_prefix}_cluster_counts.xvg"

   echo "@    title       \"Number of clusters over time\"" >${output_prefix}_nclusters.xvg
   echo "@    xaxis label \"Frame\""                       >>${output_prefix}_nclusters.xvg
   echo "@    yaxis label \"Number of clusters found\""    >>${output_prefix}_nclusters.xvg

   echo "@    title       \"Population of clusters\""       >${output_prefix}_cluster_counts.xvg
   echo "@    xaxis label \"Cluster\""                     >>${output_prefix}_cluster_counts.xvg
   echo "@    yaxis label \"Population\""                  >>${output_prefix}_cluster_counts.xvg

   for i in `seq 0 $((${#CLIDs[@]}-1))`; do
      CLID=${CLIDs[$i]}
      if ! [[ ${Unique_CLIDs[@]} =~ $CLID"," ]]; then
         Unique_CLIDs+=("$CLID,")
      fi

      echo -e "$i\t${#Unique_CLIDs[@]}" >> ${output_prefix}_nclusters.xvg

      CLID_count[$CLID]=$((${CLID_count[$CLID]}+1))
   done

   for i in `seq 1 $((${#CLID_count[@]}))`; do
      echo -e "$i\t${CLID_count[$i]}" >> ${output_prefix}_cluster_counts.xvg
   done

fi
