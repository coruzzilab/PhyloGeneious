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
#PBS -N runmclfam

echo hello for runmclfam
cd $OID_USER_DIR
OID_DATA=$OID_USER_DIR/data
MYSTART=.$PBS_JOBID.start
date +%s >$MYSTART
echo runmclfam.pl >>$MYSTART

$OID_HOME/bin/runmclfam.pl "$arg1" #$ENV_WRAPPER 
