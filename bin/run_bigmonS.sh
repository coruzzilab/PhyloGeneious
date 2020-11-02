#!/bin/bash
#
# SLURM job script for pooltnt
#
#SBATCH -o log/job/Bigmon.%J.out


echo hello for bigmon
OID_DATA=$OID_USER_DIR/data
MYSTART=.$SLURM_JOB_ID.start
$OID_HOME/bin/gettime.pl >$MYSTART
echo Bigmon.pl >>$MYSTART

cd $OID_USER_DIR
$OID_HOME/bin/Bigmon.pl
 

