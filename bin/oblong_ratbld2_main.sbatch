#!/bin/bash

#SBATCH --job-name=bgm_m2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=02:00:00


#this will need to accept arguments from preceding script:
#nchar, ntsk, iter, prefix, time from setup_proc.pl
#time = 00:30:00
#this will be called multiple times in parallel
#cd family directory (from Bigmon.pl)

if [ "$SLURM_ARRAY_TASK_ID" == "" ]; then exit; fi

num_batch=$1

source $OID_USER_DIR/family_list.sh
fams=${family_list["ratbld2"]}
read -ra families <<< "$fams"

for ((i=0;i<=num_batch;i++))
do
    index=$((SLURM_ARRAY_TASK_ID * (num_batch + 1) + i))
	if [ $index -lt ${#families[*]} ]; then
		bash $OID_HOME/bin/oblong_ratbld2_main.sh ${families[$index]}
	fi
done
