#!/bin/bash
#PBS -l nodes=1:ppn=1,pvmem=8000mb,walltime=05:00:00
#PBS -j oe
#PBS -V

module load ompi_threads/1.6.4/intel/13.1

if [ $# -lt 6 ]; then
    echo "Usage: $0 [PDB file] [dir] [IDX_file] [clust_group] [clust_cut] [identical]"
    exit
fi

################################
#  Settings from command-line  #
################################

PDB_file=$1
dir=$2
IDX_file=$3
clust_group=$4
clust_cut=$5
identical=$6

mkdir -p $dir

################
#  Clustering  #
################

time g_cluster -s ${PDB_file} -f ${PDB_file} -method gromos -n ${IDX_file} -cl $dir/clusters.pdb -cutoff ${clust_cut} -clid $dir/clust-id -o $dir/rmsd-clust.xpm -g $dir/clust.log -dist $dir/rmsd-dist -nidentical $identical >/dev/null <<< $clust_group'
4
'

###########################
#  Cluster count vs time  #
###########################

CLIDs=(`tail +19 ${dir}/clust-id.xvg | sed 's/^ *[0-9]\+ \+//'`)
Unique_CLIDs=()
N_uniq_CLIDs=0
CLID_count=()

echo "Calculating number of clusters over time in file ${dir}/nclusters.agr"
echo "Calculating population of clusters in file ${dir}/cluster_counts.agr"

echo "@    title       \"Number of clusters over time\"" >${dir}/nclusters.agr
echo "@    xaxis label \"Frame\""                       >>${dir}/nclusters.agr
echo "@    yaxis label \"Number of clusters found\""    >>${dir}/nclusters.agr

echo "@    title       \"Population of clusters\""       >${dir}/cluster_counts.agr
echo "@    xaxis label \"Cluster\""                     >>${dir}/cluster_counts.agr
echo "@    yaxis label \"Population\""                  >>${dir}/cluster_counts.agr

for i in `seq 0 $((${#CLIDs[@]}-1))`; do
    CLID=${CLIDs[$i]}

    echo -e "$i\t${#CLID_count[*]}" >> ${dir}/nclusters.agr

    CLID_count[$CLID]=$((${CLID_count[$CLID]}+1))
done

for i in `seq 1 $((${#CLID_count[@]}))`; do
    echo -e "$i\t${CLID_count[$i]}" >> ${dir}/cluster_counts.agr
done
