#!/usr/bin/python

import math, numpy, sys

kb=1.9858775e-3

class Replica:

    def __init__(self,repFile):


        """
        Initialise data for PMF window
        """
        
        ## store file name
        self.fname=repFile

        ## open file containing probability data
        dataF=open(repFile,'r')
        lines=dataF.readlines()

        ## get no of bins
        self.nbin=len(lines)-1

        ## read in k and x0
        data=lines[0].split()
        self.temp=float(data[1])
        self.beta=1.0/(kb*self.temp)

        # initialise energy and probablilty arrays
        self.eng=numpy.zeros((self.nbin),dtype=numpy.float128)
        self.hist=numpy.zeros((self.nbin),dtype=numpy.float128)
        s=0.0
        for i,line in enumerate(lines[1:]):

            [e,p]=[float(d) for d in line.split()]

            self.eng[i]=e
            self.hist[i]=p


        ## set initial partition function
        self.z=1.0
        self.zOld=1.0


def whamIter(replicas):

    nRep=len(replicas)

    ## loop over replicas and save current parition functions
    for replica in replicas:
        replica.zOld=replica.z
        replica.z=0.0
    
    for s in range(replicas[0].nbin):
        sum_i=0.0
        eng=replicas[0].eng[s]
        for i in range(nRep):
            sum_i+=replicas[i].hist[s]
        for k in range(nRep):
            sum_j=0.0
            for j in range(nRep):
                sum_j+=(1.0/replicas[j].zOld)*numpy.exp((replicas[k].beta-replicas[j].beta)*eng)
            replicas[k].z+=sum_i/sum_j
        
    


if __name__=='__main__':

    ## get arguments from command line

    ## first check for number of agruments
    if len(sys.argv) < 4:

        print "Error: incorrect no of arguments"
        sys.exit()
        
    ## get tolerance and max iterations off end of argument list
    tol=float(sys.argv.pop())
    maxit=int(sys.argv.pop())
    outFile=sys.argv.pop()
        
    ## get no of input files
    nfile=len(sys.argv[1:])

    ## initialise list of pmf data
    replicas=[]

    ## read in data for each window
    for f in sys.argv[1:]:

        print "## reading data from file ",f
        replicas.append(Replica(f))

    for i in range(maxit):
        whamIter(replicas)
        print i+1,
        delta=0.0
        for rep in replicas:
            print rep.z,
            delta+=((rep.z-rep.zOld)/rep.z)**2
        print delta
        if delta<tol:
            break
        
    else:
        print "No of iterations exceeded"

    ouf=open(outFile,'w')
    print "## final partition functions"
    for i,replica in enumerate(replicas):
        print "## replica %6d partition function %12.6g " % (i+1,replica.z)

    print >> ouf, "## final partition functions"
    for i,replica in enumerate(replicas):
        print >> ouf, "## replica %6d temperature %12.6f partition function %12.6g " % (i+1,replica.temp,replica.z)


    ## create final histogram
    nbin=replicas[0].nbin
    histFinal=numpy.zeros((replicas[0].nbin),dtype=numpy.float128)
    for s in range(nbin):
        sum_i=0.0
        sum_j=0.0
        e=replicas[0].eng[s]
        for i in range(len(replicas)):
            sum_i+=replicas[i].hist[s]
        for j in range(len(replicas)):
            sum_j+=(1.0/replicas[j].z)*numpy.exp(-replicas[j].beta*e)
        print s,sum_i,sum_j
        histFinal[s]=sum_i/sum_j

    norm=sum(histFinal)
    for i in range(nbin):
        small=1.0e-100
        if histFinal[i]<=small:
            continue
        print >> ouf, "%12.6f %12.6g " % (replicas[0].eng[i],numpy.log(histFinal[i]))
