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
##PBS -o /scratch/cmz209/orthotnt/oidTest8/log/job
#PBS -q cgsb-s
#PBS -N runmclfam
##PBS -M cmz209@nyu.edu
##PBS -m abe


echo hello for runmclfam
cd $OID_USER_DIR
OID_DATA=$OID_USER_DIR/data
MYSTART=.$PBS_JOBID.start
$OID_HOME/bin/gettime.pl >$MYSTART
echo runmclfam.pl >>$MYSTART

$OID_HOME/bin/runmclfam.pl "$arg1"
