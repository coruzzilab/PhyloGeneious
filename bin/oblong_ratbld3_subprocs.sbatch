#!/bin/bash

#SBATCH --job-name=bgm_s3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=12:00:00


#this will need to accept arguments from preceding script:
#nchar, ntsk, iter, prefix, time from setup_proc.pl
#time = 00:30:00
#this will be called multiple times in parallel
#cd family directory (from Bigmon.pl)

if [ "$SLURM_ARRAY_TASK_ID" == "" ]; then exit; fi

num_batch=$1
subprocs_ntimes=10

source $OID_USER_DIR/family_list.sh
fams=${family_list["ratbld3"]}
read -ra families <<< "$fams"

num_fams=${#families[@]}

# ntimes=10

for ((i=0;i<=num_batch;i++))
do
    index=$((SLURM_ARRAY_TASK_ID * (num_batch + 1) + i))
    if [ $index -lt $((num_fams * subprocs_ntimes)) ]; then
        bash $OID_HOME/bin/oblong_ratbld3_subprocs.sh ${families[$((index % $num_fams))]} $((index / $num_fams))
    fi
done
