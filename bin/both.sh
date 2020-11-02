#!/bin/bash
#SBATCH -p parallel
#SBATCH --mem 22GB
#SBATCH --ntasks=12
#SBATCH --natasks-per-core=1
##SBATCH  -c 4
#SBATCH  -cpus-per-task=4
#SBATCH --time=02:30:00
#SBATCH -o job.%J.out

# **************
echo hello
#srun -p serial --mem 22GB -n 1 -c 8 -o job.%J.out $OID_HOME/bin/run_just_bl.sh&
export MY_SHELL=$0
sbatch -c 4 --mem 22GB  -t 04:00:00 -o toplog/%J.out $OID_HOME/bin/run_both.sh

