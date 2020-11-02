#!/bin/bash
#
# SLURM job script for pooltnt
#
#SBATCH -o log/job/pooltnt.%J.out


cd $OID_USER_DIR
MYSTART=.$SLURM_JOB_ID.start
$OID_HOME/bin/gettime.pl >$MYSTART
echo pooltnt.pl $1 $2 $3 >>$MYSTART
date
time
$OID_HOME/bin/pooltnt.pl $1 $2 $3

date
time
