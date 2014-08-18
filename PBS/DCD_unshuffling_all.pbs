#!/bin/bash
#PBS -l nodes=1:ppn=1,pvmem=2000mb,walltime=72:00:00
#PBS -j oe
#PBS -V

##############
#  Settings  #
##############

first_dump_step=0
freq_dump_steps=10000

####################################################
#  Sanity check: Is this job part of a job array?  #
####################################################

number_regex='^[0-9]+$'
if ! [[ ${PBS_ARRAYID} =~ $number_regex ]]; then
   echo "Error: This isn't an array"
   exit 1
fi

cd $PBS_O_WORKDIR

###############################################################
#  Unshuffle, i.e. make trajectories of constant temperature  #
###############################################################

#Run as: ~/scripts/post_process_md/REMD_dcd_unshuffler.py [filestring] [F] [N] [T] {B}
#[filestring] If the files are called traj_0.dcd .. traj_N.dcd, filestring is 'traj_'.
#[F] First dump step.
#[N] Frequency of dump steps.
#[T] Which temperature to demultiplex
#{B} Which DCD frame to start with. Typically set this to 1 more than the number of frames so far.  = 1 by default.
   #Run in a directory containing the files to be unshuffled and the log.lammps corresponding to them.

$MDdir/REMD/dcd_unshuffler.py '' ${first_dump_step} ${freq_dump_steps} ${PBS_ARRAYID}
rename '-1' '' *.dcd

##########################################
#  Unwrap molecules to make pbcs better  #
##########################################

module unload vmd-1.9.1
module load vmd-1.8.6

vmd T${PBS_ARRAYID}.dcd -psf ../*.psf -dispdev text -e $MDdir/unwrap_script.tcl -args ${PBS_ARRAYID}
