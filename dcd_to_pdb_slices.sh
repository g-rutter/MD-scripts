#!/bin/sh

if ! [ $# -eq 4 ]; then
    echo "Usage: $0 [DCD_file] [PSF_file] [stride] [percent_slice_size]"
    exit
fi

DCD=$1
PSF=$2
stride=$3
percent_slice_size=$4
threads=2

start_percents=(`seq 0 $percent_slice_size 99`)

end_percents=(`seq $((0+percent_slice_size)) $percent_slice_size 99` 100)

echo ${start_percents[*]}
echo ${end_percents[*]}

echo ${#start_percents[*]}
echo ${#end_percents[*]}

for i in `seq 0 $(( ${#end_percents[*]} - 1 ))`; do
    start=${start_percents[$i]}
    end=${end_percents[$i]}

    echo $start $end
    $MDdir/dcd_to_pdb.sh $DCD $PSF $stride $start $end
done
