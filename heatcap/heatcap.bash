#!/bin/bash
#Try: for i in {0..44}; do   num=`sed '738q;d' log.lammps.$i|sed 's/ \([0-9]*\).*/\1/'`; if [ $num != 1001000 ]; then sed -i '1s/^/\n/' log.lammps.$i; fi; done
#Use this to pad log.lammps files if they cause an inconsistency error with thermo_table_maker

#Read-in args
args=("$@")
n_args=${#args[*]}
n_temps=$(($n_args-2))
el=${args[0]}
eh=${args[1]}

if [ $n_args = 2 ]; then
   echo E range: $el $eh
else
   echo "Please pass 2 command line arguments; the low and high energy range bounds."
   exit
fi

#In-file settings
script_path="$HOME/scripts/post_process_md"
log_prefix="log.lammps."
temp_range=(270 400)
outfile_dir="heatcapdata"

lmp_inp_search="../*.in"
lmp_inp=(`ls $lmp_inp_search`)

if ! [ ${#lmp_inp[*]} = 1 ]; then
   echo WARNING, found not one lmp input files: ${lmp_inp[*]}
   echo Using $lmp_inp
fi

temps=(`ack -i '^variable t world' $lmp_inp | sed 's/^variable t world *//'`)
n_temps=${#temps[*]}
echo ${temps[*]}

############################
# Prepare for file outputs #
############################

mkdir -p $outfile_dir

if ls $outfile_dir/* >/dev/null 2>&1; then
   echo "Delete your old files in $outfile_dir if your input files have changed."
fi

#####################
# Make TotEng table #
#####################

if ! [ -f $outfile_dir/"TotEngs" ]; then
   if [ -f ${log_prefix}0 ]; then
      first_line=$((`ack '^Step.*TotEng' log.lammps.0 --heading -H | tail -1 | sed 's/:.*$//'`+1))
      columns=(`ack '^Step.*TotEng' log.lammps.0 --heading -H | tail -1`)
      for (( i = 0; i < ${#columns[*]}; i++ )); do 
         if [ "${columns[i]}" = "TotEng" ]; then
            CV_column=$(($i+1))
         fi
      done
   else
      echo "Can't find "$log_prefix"0"
      exit
   fi

   ${script_path}/REMD_thermo_table_maker.py $log_prefix $first_line 1 $CV_column
   mv CVs $outfile_dir/TotEngs
fi

#####################
# Make E histograms #
#####################

if ! [ -f $outfile_dir/E_bounds ]; then
   echo -e "999999\n-999999" >${outfile_dir}/E_bounds
fi
echo -ne "" >${outfile_dir}/Emins
echo -ne "" >${outfile_dir}/Emaxs
echo -ne "" >${outfile_dir}/Eavgs

echo Histogramming temperature data

for i in $(seq 1 $n_temps); do
   T=${temps[$(($i-1))]}
   if ! [ -f $outfile_dir/T${T}.hist ]; then
      E_hist_out="`python ${script_path}/heatcap/makeEHist.py $outfile_dir/TotEngs $outfile_dir/T${T}.hist $el $eh 1500 ${T} $i`"

      echo ""
      echo "$E_hist_out"
      echo ""
      #Keep track of min and max E values.
      new_min=`echo "$E_hist_out" | ack "min" | sed 's/^min \+= \+//'`
      new_max=`echo "$E_hist_out" | ack "max" | sed 's/^max \+= \+//'`
      new_avg=`echo "$E_hist_out" | ack "average" | sed 's/^average \+= \+//'`
      old_min=`head -1 ${outfile_dir}/E_bounds 2>/dev/null `
      old_max=`tail -1 ${outfile_dir}/E_bounds 2>/dev/null `

      if [ `echo "$new_min < $old_min" | bc -l` -eq 1 ]; then
         old_min=$new_min
      fi
      if [ `echo "$new_max > $old_max" | bc -l` -eq 1 ]; then
         old_max=$new_max
      fi

      echo -e "$old_min\n$old_max" >${outfile_dir}/E_bounds
      echo -e "$T $new_min" >>${outfile_dir}/Emins
      echo -e "$T $new_max" >>${outfile_dir}/Emaxs
      echo -e "$T $new_avg" >>${outfile_dir}/Eavgs

   fi
done

echo Minimum energy $old_min
echo Maximum energy $old_max

if ! [ -f $outfile_dir/log_rho.dat ]; then
   python ${script_path}/heatcap/ptWhamHist.py $outfile_dir/T*.hist $outfile_dir/log_rho.dat 100000 1.0e-8
fi

python ${script_path}/heatcap/cv.py $outfile_dir/log_rho.dat heatcap ${temp_range[*]} 1
