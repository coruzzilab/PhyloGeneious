#!/bin/bash
#
# PBS job script for orthologid
#

#PBS -V
#PBS -S /bin/bash
#PBS -j oe
#PBS -o log/job
#PBS -l mem=12GB
#PBS -l nodes=1:ppn=1,walltime=48:00:00
#PBS -N runmcx

echo hello for runmcx
cd $OID_USER_DIR
OID_DATA=$OID_USER_DIR/data
MYSTART=.$PBS_JOBID.start
date +%s >$MYSTART
echo runmcx.pl >>$MYSTART

$OID_HOME/bin/runmcx.pl #$ENV_WRAPPER 
