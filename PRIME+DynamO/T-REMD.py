#!/usr/bin/env python
# encoding: utf-8

from argparse import ArgumentParser
import bz2
import math
import numpy
import os
import random
import re
import shutil
import subprocess
import time

##########
#  Help  #
##########

#Important variables:
#T_index[i_replica] = Temperature_index of replica i

#######################
#  In-script options  #
#######################

valid_modes = ["_20-like", "_2001", "_2001-extendedHB"]
dynasuffix="_2001-extendedHB"

main_files = 'run_files/'


########################
#  Command-line input  #
########################

#Set up argument parser
parser = ArgumentParser()
parser.add_argument("infiles", nargs='+', help="Files to start the simulation from, each held at a different temperature.")
parser.add_argument('-s', type=int, help="Swap frequency; time units.")
parser.add_argument('-r', type=int, help="Total run time; time units.")
parser.add_argument('-m', type=str, help="Dynarun/dynamod mode. Default="+dynasuffix+". Options: "+", ".join(valid_modes))

#Interpret input
args  = parser.parse_args()

infiles=args.infiles
swap_time_root=args.s
run_time_root=args.r

if args.m is not None:
    if args.m in valid_modes:
        dynasuffix=args.m
    else:
        print "-m needs to be one of:", ", ".join(valid_modes)
        exit()

dynarun="dynarun"+dynasuffix
dynamod="dynamod"+dynasuffix

############################
#  Read files, get set-up  #
############################

Ts=[]
run_files=[]
temp_expr = re.compile("Temperature=\"([0-9\.]+)\"")
n_atoms_expr = re.compile("<Pt ID=")

N_atoms = None

if not os.path.exists(main_files):
    os.mkdir(main_files)

for infile in infiles:
    infile_ptr      = bz2.BZ2File(infile)
    infile_contents = infile_ptr.read()
    Ts.append( float(temp_expr.search(infile_contents).group( 1 )) )

    run_files.append(main_files+str(len(run_files))+".run.xml.bz2")
    shutil.copy(infile, run_files[-1])

    N_atoms_temp = len ( n_atoms_expr.findall(infile_contents) )
    if N_atoms is None:
        N_atoms = N_atoms_temp
    else:
        try:
            assert (N_atoms == N_atoms_temp)
        except AssertionError:
            print "Files don't contain same atom counts. Make sure the files are the same system."
            exit()

N_replicas = len(Ts)
N_swaps    = int( math.ceil( float(run_time_root)/swap_time_root ) )
LowT       = min(Ts)
swap_times = []

with open(main_files + 'log.DynamO', 'w') as REMD_log:
    REMD_log.write("Step "+" ".join(['T'+str(i) for i in range(N_replicas)]) + '\n')
    REMD_log.write("0 "+" ".join([str(i) for i in range(N_replicas)]) + '\n')

for T in Ts:
    swap_times.append(str(math.sqrt(LowT/T)*swap_time_root))

for j_temperature in range(N_replicas):
    T_str=str(j_temperature)
    if not os.path.exists("T"+T_str):
        os.mkdir("T"+T_str)

T_index = numpy.array( range(N_replicas) )
swap_generator = list( range(N_replicas) )

#Swapping success trackers
swap_tries = numpy.zeros([N_replicas,N_replicas], dtype=int)
swap_successes = numpy.zeros([N_replicas,N_replicas], dtype=int)

adjacent_pairs = tuple(tuple([i,i+1]) for i in range(0,N_replicas-1,2))

print "Temperatures:", Ts

#########
#  Run  #
#########

def service_finished_runs(processes):

    N_finished_processes = 0
    Us=numpy.empty(N_replicas, dtype=float)
    ETA_expr = re.compile("ETA 0s, .*, (T .*, U ([-0-9\.]+))")
    HBond_expr = re.compile("Interaction\[Backbone.*")

    filename="/Snapshot."+str(i_swap)+".xml.bz2"

    while N_finished_processes != N_replicas:

        for j_replica, process in enumerate(processes):
            #Pass if already handled
            if process is None:
                continue

            #Handle newly finished
            if process.poll() is not None:
                HBond_matches = HBond_expr.findall( process.communicate()[0] )
                processes[j_replica]=None
                N_finished_processes+=1

                shutil.copy(run_files[j_replica], "T"+str(T_index[j_replica]) + filename)

                match = ETA_expr.search(
                                        subprocess.Popen(
                                            [ dynarun, "-c", "1", "-o", "/dev/null", run_files[j_replica] ],
                                            stdout=subprocess.PIPE
                                        ).communicate()[0]
                                    )

                with open("T"+str(T_index[j_replica])+"/ETA", 'a') as ETAfile:
                    ETAfile.write( match.group(1)+'\n' )
                Us[j_replica]=float ( match.group(2) ) * N_atoms

                with open("T"+str(T_index[j_replica])+"/HB", 'a') as HBfile:
                    HBfile.write( '\n'.join(HBond_matches) )
                    HBfile.write( '\n\n\n' )

        time.sleep(.03)

    return Us

def bool_swap(T_i,T_j,U_i,U_j):

    #print "U_"+str(T_i)+" =", str(U_i)
    #print "U_"+str(T_j)+" =", str(U_j)

    dU=U_i-U_j
    dBeta=(1/T_i)-(1/T_j)
    p=math.exp(dU*dBeta)

    #print "delta beta =", str(dBeta)
    #print "delta U    =", str(dU)
    #print "p          =", str(p)

    if p>1:
        return True

    rand = random.random()
    #print "Random number [0,1]:", rand

    if rand < p:
        return True

    return False

for i_swap in range(N_swaps):

    #Submit subruns
    processes = list()
    for j_replica in range(N_replicas):
        processes.append(subprocess.Popen(
            [dynarun, run_files[j_replica], '-f', swap_times[j_replica], '-o', run_files[j_replica], '-s', '1'],
            stdout=subprocess.PIPE
        ))

    step=str((1+i_swap)*swap_time_root)+" "
    random.shuffle(swap_generator)

    #Let processes terminate and get final potentials indexed by replica
    Us = service_finished_runs(processes)
    #print ""
    #print "{0:.2f}% complete.".format( 100*(1.0+i_swap)/N_swaps )

    #Produce list of swaps to attempt
    for a, b in adjacent_pairs:
        i_replica = swap_generator[a]
        j_replica = swap_generator[b]

        T_i=Ts[T_index[i_replica]]
        T_j=Ts[T_index[j_replica]]

        swap_tries[T_index[i_replica]][T_index[j_replica]] += 1
        swap_tries[T_index[j_replica]][T_index[i_replica]] += 1

        if bool_swap(T_i, T_j, Us[i_replica], Us[j_replica]):

            swap_successes[T_index[i_replica]][T_index[j_replica]] += 1
            swap_successes[T_index[j_replica]][T_index[i_replica]] += 1

            #print "swap! {0}({1}) {2}({3})".format( i_replica, T_i, j_replica, T_j )
            T_index[i_replica], T_index[j_replica] = T_index[j_replica], T_index[i_replica]

            subprocess.Popen(
                [ dynamod, "-r", str(T_j), "-T", str(T_j), run_files[i_replica], '-o', run_files[i_replica] ],
            stdout=subprocess.PIPE).communicate()

            subprocess.Popen(
                [ dynamod, "-r", str(T_i), "-T", str(T_i), run_files[j_replica], '-o', run_files[j_replica] ],
            stdout=subprocess.PIPE).communicate()

    #for j_replica in range(N_replicas):
        #print "Replica", j_replica, "should now be at temperature", Ts[T_index[j_replica]]
    #print ""

    with open(main_files+'log.DynamO', 'a') as REMD_log:
        REMD_log.write(step+" ".join([str(numpy.where(T_index == i)[0][0]) for i in range(N_replicas)])+'\n' )

#############################
#  Output some REMD stats.  #
#############################

swap_ratio = numpy.nan_to_num( numpy.true_divide( swap_successes, swap_tries) )

percent_formatstr=" ".join(["{"+str(i)+":5.1f}" for i in range(N_replicas)])
temp_formatstr=" ".join(["{"+str(i)+":5}" for i in range(N_replicas)])
print "{0:5}".format("Temp"), temp_formatstr.format(*Ts)

for i in range(N_replicas):
    print "{0:5}".format(Ts[i]),
    print percent_formatstr.format(*100*swap_ratio[i])
