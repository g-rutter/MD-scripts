#!/bin/sh

##############
#  Settings  #
##############

unite_dir="run_all"

############
#  Set-up  #
############

if [ -d $unite_dir ]; then
    echo "Warning: $unite_dir exists."
else
    mkdir $unite_dir
fi

N_replicas=`find "run_1" -maxdepth 1 -name '[0-9]*' | wc -l`
run_dirs=(run_[0-9]*)
N_runs=${#run_dirs[*]}

echo "Found $N_runs run_dirs, $N_replicas replicas."

j_replicas=($*)
if [ -z "$j_replicas" ]; then
    echo Sorting all replicas
    j_replicas=(`seq 0 $(( N_replicas - 1 ))`)
else
    echo Sorting replicas: ${j_replicas[*]}
fi

#############################
#  Copy and join snapshots  #
#############################

for j_replica in ${j_replicas[*]}; do

    mkdir -p $unite_dir/$j_replica
    snap_tally=0
    echo -ne "Replica $j_replica: "

    for i_run in `seq 1 $N_runs`; do
        echo -ne "$i_run "
        snaps=(`find "run_$i_run/$j_replica" -name 'Snapshot.ID*'`)
        N_snaps_here=${#snaps[*]}
        for k_snap in `seq 0 $(( N_snaps_here - 1 ))`; do

            l_snap=$((snap_tally+k_snap))
            this_snap=${snaps[$k_snap]}

            cp $this_snap $unite_dir/$j_replica/Snapshot.ID$j_replica.$l_snap.xml.bz2
        done

        snap_tally=$((snap_tally + N_snaps_here))
    done

    echo ""
done
