#!/bin/bash

#SBATCH --job-name=run_mafft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=00:30:00


if [ "$SLURM_ARRAY_TASK_ID" == "" ]; then exit; fi

readarray -t families < $OID_USER_DIR/FAMILY

thread=4
num_batch=$1

for ((i=0;i<=num_batch;i++))
do
    index=$((SLURM_ARRAY_TASK_ID * (num_batch + 1) + i))
    if [[ $index -lt ${#families[*]} ]]; then
        bash $OID_HOME/bin/run_mafftS.sh ${families[$index]} $thread
    fi
done
