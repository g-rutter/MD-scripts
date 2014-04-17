#!/usr/bin/python

import numpy, sys

kb=1.9858775e-3

inf=open(sys.argv[1],'r')
ouf=open(sys.argv[2],'w')

temp1=float(sys.argv[3])
temp2=float(sys.argv[4])
dtemp=float(sys.argv[5])

lines=inf.readlines()

eng=[]
logrho=[]
for line in lines:
    if line[0]=='#':
        continue
    data=line.split()
    eng.append(float(data[0]))
    logrho.append(float(data[1]))

eng=numpy.array(eng,dtype=numpy.float128)
logrho=numpy.array(logrho,dtype=numpy.float128)

temp=temp1
while 1:

    beta=1.0/(kb*temp)
    ave_e=0.0
    ave_e2=0.0

    logp=logrho-beta*eng
    p=numpy.exp(logp)
    sum_p=sum(p)
    p/=sum_p

    ave_e=sum(eng*p)
    ave_e2=sum(eng**2*p)

    cv=ave_e2-ave_e**2

    print >> ouf, temp,cv
    print temp,cv

    temp+=dtemp
    if temp>temp2:
        break
