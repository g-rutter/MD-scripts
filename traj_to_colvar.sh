#!/bin/bash

DCDs=(`find . -name "T*.dcd" | sort -n -t T -k 2`)
N_trajs=${#DCDs[*]}

PLUMED=(`find .. -maxdepth 2 -name *.inp`)
PDB=(`find .. -maxdepth 2 -name *.pdb`)
PSF=(`find .. -maxdepth 2 -name *.psf`)

echo Trajectory files: ${DCDs[*]}
echo PLUMED: ${PLUMED[*]}
echo PDB: ${PDB[*]}

if [ ${#PLUMED[*]} -gt 1 ]; then
   echo "WARNING: found >1 PLUMED .inp files: ${PLUMED[*]}"
   echo Using $PLUMED
elif [ ${#PLUMED[*]} -eq 0 ]; then
   echo ERROR: Found 0 PLUMED files.
   exit
fi

if [ ${#PDB[*]} -gt 1 ]; then
   echo "WARNING: found >1 PDB files: ${PDB[*]}"
   echo Using $PDB
elif [ ${#PDB[*]} -eq 0 ]; then
   echo ERROR: Found 0 PDB files.
   PDB="1.pdb"
   echo catdcd -last 1 -o $PDB -otype pdb -stype  -s $PSF -dcd ${DCDs[0]}
   catdcd -last 1 -o $PDB -otype pdb -stype psf -s $PSF -dcd ${DCDs[0]}
fi

#Ask user to verify settings
echo "Proceed? [y/n]"
read A
if ! [ $A = "y" -o $A = "Y" ]; then
   exit
fi

echo ""
echo -n "Silently running gpta.x on temperature index: "

#Print colvars to COLVAR-T$i files.
for (( T = 0; T < $N_trajs; T++ )); do
   command="gpta.x -ipdb $PDB -bidcd ${DCDs[$T]} -plumed $PLUMED"
   echo -n "${T} "
   $command >/dev/null

   awk -F " *" '{print $3}' COLVAR >COLVAR.T$T
   COLVARs[$T]=COLVAR.T${T}
done

echo ""
echo "Combining files."

paste ${COLVARs[*]} | tail --lines=+2 >"COLVAR_all"

rm -f PLUMED.OUT
rm -f COLVAR.T* COLVAR COLVAR.old
