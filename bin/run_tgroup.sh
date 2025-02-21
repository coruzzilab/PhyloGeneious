#!/bin/sh
#
# PBS job script for orthologid
#

#PBS -S /bin/bash
#PBS -j oe
#PBS -o log/job
#PBS -l mem=12GB
#PBS -N pooltnt
#PBS -V

cd $OID_USER_DIR
date
time
$OID_HOME/bin/pooltnt.pl "$arg1" "$arg2" "$arg3"  #$ENV_WRAPPER 

date
time
