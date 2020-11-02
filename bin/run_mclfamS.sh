#!/bin/bash
#
# SLURM job script for runmclfam.pl
#
#SBATCH -o log/job/mclfam.%J.out


echo hello for runmclfam
OID_DATA=$OID_USER_DIR/data
MYSTART=.$SLURM_JOB_ID.start
$OID_HOME/bin/gettime.pl >$MYSTART
echo runmclfam.pl >>$MYSTART

cd $OID_USER_DIR
$OID_HOME/bin/runmclfam.pl $1
 

