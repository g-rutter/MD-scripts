#!/usr/bin/env python
# encoding: utf-8

from numpy import array, empty
from sys import argv

######################
#  In-fn settings  #
######################

units = 'real'
kB    = 0.001987 # Boltzmann constant in kcal/mol/K

######################
#  Read-in settings  #
######################

if len(argv) < 3:
    print "Run as:", argv[0], "[Temps fn] [TotEngs fn]"
    print "[Temps fn] file with first line containing space-delimited list of temperatures"
    print '''[TotEngs fn] file from $MDdir/REMD/thermo_table_maker.py giving space-delimited
               ordered list of total energy in each temperature, each line being a snapshot.'''
    exit()

T_fn = argv[1]
totengs_fn = argv[2]

####################
#  Validate fns  #
####################

with open(T_fn) as T_file:
    Ts = [ float(T_str) for T_str in T_file.readline().split() ]
    N_Ts = len(Ts)
    Ts = array( Ts )

with open(totengs_fn) as TE_file:
    try:
        N_Ts_check = len( TE_file.readline().split() )
        assert ( N_Ts == N_Ts_check )
    except AssertionError:
        print "According to temps file there are", N_Ts, "temperatures."
        print "But the first line in TotEngs file has", N_Ts_check, "entries."
        exit()

    N_snapshots = 1 + sum(1 for _ in TE_file)

print "Detected", N_Ts, "temperatures and", N_snapshots, "snapshots."

######################
#  Read in Tot Engs  #
######################

TotEngs = empty([N_snapshots,N_Ts], dtype=float)

with open(totengs_fn) as TE_file:
    for i_snap, snapshot_str in enumerate(TE_file):
        TotEngs[i_snap] = snapshot_str.split(' ')

################################
#  Heat cap method 1: variance #
################################

TotEng_variances = TotEngs.var(axis=0)
Cv_1 = TotEng_variances/(kB*Ts**2)
print "Variance method:", Cv_1

#################################################
#  Heat cap method 2: gradient of average TotE  #
#################################################

TotEng_avgs = TotEngs.mean(axis=0)

method2_Ts = empty([N_Ts-1], dtype=float)
Cv_2 = empty([N_Ts-1], dtype=float)

for i in range(0, N_Ts-1):
    method2_Ts[i] = Ts[i: i+2].mean()
    Cv_2[i]       = (TotEng_avgs[i+1] - TotEng_avgs[i])/(Ts[i+1]-Ts[i])

print "Gradient method:", Cv_2

#############
#  xmgrace  #
#############

with open('variance.agr', 'w') as varfile:
    for i_temp in range(N_Ts):
        varfile.write( str(Ts[i_temp]) + ' ' + str(Cv_1[i_temp]) + '\n' )

with open('gradient.agr', 'w') as gradfile:
    for i_temp in range(N_Ts-1):
        gradfile.write( str(method2_Ts[i_temp]) + ' ' + str(Cv_2[i_temp]) + '\n' )
