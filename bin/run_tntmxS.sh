#!/bin/bash
#
# SLURM job script for pooltnt
#
#SBATCH -o log/job/tntmx.%J.out


echo hello for $1
OID_DATA=$OID_USER_DIR/data
MYSTART=.$SLURM_JOB_ID.start
$OID_HOME/bin/gettime.pl >$MYSTART
echo runtntmx.pl $1 $2 $3 $4 >>$MYSTART
fam=$1
cd $OID_DATA/$fam
echo `pwd`
echo $2 $3 $4
$OID_HOME/bin/runtntmx.pl $2 $3 $4

 

