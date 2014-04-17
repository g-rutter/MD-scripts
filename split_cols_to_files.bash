#!/bin/bash

FILE=$1
firstline=(`head -1 $FILE`)
n_cols=${#firstline[*]}

echo "It looks like there are "$n_cols" columns."

echo rm -I -v $FILE.*
rm -I -v $FILE.*
echo "Splitting to separate files."

for j_col in `seq 1 $n_cols`; do

   cut -f $j_col -d ' ' $FILE >$FILE.$j_col

done
