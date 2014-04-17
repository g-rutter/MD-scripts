## script to contruct energy histogram
#!/usr/bin/python

import gzip, numpy, sys

def makeHist(xlo,xhi,nx,data):

    dx=(xhi-xlo)/nx
    hist=numpy.zeros((nx))
    for d in data:
        if d<xlo or d>xhi:continue
        dd=d-xlo
        bin=int(numpy.floor(dd/dx))
        hist[bin]+=1.0

    return hist

if len(sys.argv) != 8:
   print sys.argv
   print "Run as: ", sys.argv[0], "[datafile] [histfile] [Elo] [Ehi] [n Bins] [temp] [E_tot column]"
   print "[datafile] contains one sample per line, with the total energy being the [E_tot column]th entry"
   print "[histfile] output filename"
   print "[Elo] [Ehi] bounds for energy histogram"
   print "[n bins] number of bins"
   print "[temp] temperature data corresponds to"
   exit()
dataFile=sys.argv[1]
histFile=sys.argv[2]
elo=float(sys.argv[3])
ehi=float(sys.argv[4])
nBin=int(sys.argv[5])
temp=float(sys.argv[6])
e_tot_column=int(sys.argv[7])

eSeries=[]

if dataFile[-3:]=='.gz':
    df=gzip.open(dataFile,'r')
else:
    df=open(dataFile,'r')

eAve=0.0
eSqr=0.0
eMin=1.0e10
eMax=-1.0e10
cnt=0

## read in data
while 1:
    line=df.readline().strip()
    if line=='':
        break

    if line[0]=='#':
        continue

    data=line.split()
    e=float(data[e_tot_column-1])
    eSeries.append(e)

    cnt+=1
    eAve+=e
    eSqr+=e**2

    if e<eMin:
        eMin=e
    if e>eMax:
        eMax=e

eAve/=cnt
eSqr/=cnt
eStd=numpy.sqrt(eSqr-eAve**2)
print "average = %12.6f " % (eAve)
print "std dev = %12.6f " % (eStd)
print "min     = %12.6f " % (eMin)
print "max     = %12.6f " % (eMax)

## get histogram
eHist=makeHist(elo,ehi,nBin,eSeries)/cnt
hf=open(histFile,'w')
print >> hf, "#### %12.6f " % (temp)
de=(ehi-elo)/nBin
for i,h in enumerate(eHist):
    e=elo+(i+0.5)*de
    print >> hf, "%12.6f %12.6f " % (e,h)
