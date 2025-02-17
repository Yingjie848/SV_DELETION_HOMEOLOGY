#!/usr/bin/bash
# Description: submit job to LSF cluster
# Author: Yingjie Zhu

# monitor progress by examing file generated or not
function monitor_progress {

  echo "monitoring progress"
  in_file=$1
  duration_hours=72
  interval_seconds=30
  n=$(($duration_hours*3600/$interval_seconds))
  for (( i=1; i<=$n; i++ ))
  do
    echo -n "Monitoring task $i "; date
    if [[ -f $in_file ]]
    then
      echo "Job finished!"
      break
    fi
    sleep $interval_seconds
  done

}


walltime=$1
memory=$2
threads=$3
err_file=$4
out_file=$5
commands=$6

echo [walltime $walltime] [memory $memory] [threads $threads] [err_file $err_file] [out_file $out_file] [command $commands]

# submit by bsub
if [ -f $err_file ];then
  rm $err_file
fi

if [ -f $out_file ];then
  rm $out_file
fi

bsub -W $walltime -R "rusage[mem=$memory]" -n $threads -e $err_file -o $out_file $commands

# monitor progress
monitor_progress $out_file
