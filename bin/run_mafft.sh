#!/bin/sh
#
# PBS job script for orthologid
#

#PBS -S /bin/bash
#PBS -j oe
#PBS -l mem=12GB
#PBS -o log/job
#PBS -N bigmon
#PBS -V

cd $OID_USER_DIR
date
time

echo hello for mafft
OID_DATA=$OID_USER_DIR/data

cd $OID_DATA/"$arg1"
$ENV_WRAPPER mafft --auto --quiet --anysymbol --thread "$arg2" FAMILY >FAMILY.aligned
cd $OID_USER_DIR
$OID_HOME/bin/makenex.pl "$arg1" #$ENV_WRAPPER 

date
time
