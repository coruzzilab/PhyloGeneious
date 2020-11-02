#!/bin/bash
#
# SLURM job script for runmcl.pl
#
#SBATCH -o log/job/runmcx.%J.out


echo hello for runmcx
OID_DATA=$OID_USER_DIR/data
MYSTART=.$SLURM_JOB_ID.start
$OID_HOME/bin/gettime.pl >$MYSTART
echo runmcx.pl >>$MYSTART

cd $OID_USER_DIR
$OID_HOME/bin/runmcx.pl
 

