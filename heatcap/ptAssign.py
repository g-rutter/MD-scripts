## Script to asssign thermodynamic properties to correct temperature from LAMMPS parallel tempering simulation
#!/usr/bin/python

import sys, gzip

nReplica=int(sys.argv[1])
logF=open(sys.argv[2],'r')
repF=[]
repT=[]
repFout=[]
for i in range(nReplica):
    if sys.argv[3+i][-3:]=='.gz':
        repF.append(gzip.open(sys.argv[3+i],'r'))
    else:
        repF.append(open(sys.argv[3+i],'r'))
for i in range(nReplica):
    repT.append(float(sys.argv[3+nReplica+i]))
for i in range(nReplica):
    if sys.argv[3+2*nReplica+i][-3:]=='.gz':
        repFout.append(gzip.open(sys.argv[3+2*nReplica+i],'w'))
    else:
        repFout.append(open(sys.argv[3+2*nReplica+i],'w'))


print repT

## read in data from master log file
lines=logF.readlines()[3:]

for rr in repF:
    rr.readline()

## get initial assignent of replicas
data=lines[0].split()
repAss=[int(d) for d in data[1:]]
for line in lines[1:]:
    data=line.split()
    tNext=int(data[0])
    while 1:
        therm=[]
        for i in range(nReplica):
            tline=repF[i].readline().strip()
            
            t=int(tline.split()[0])
            therm.append(tline)
        if t<tNext:
            for i in range(nReplica):
                print >> repFout[repAss[i]], therm[i],i
        else:
            repAss=[int(d) for d in data[1:]]
            for i in range(nReplica):
                print >> repFout[repAss[i]], therm[i],i
            break
    
