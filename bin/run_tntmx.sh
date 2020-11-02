#!/bin/bash
#
# PBS job script for orthologid
#

#PBS -V 
#PBS -S /bin/bash
#PBS -j oe
##PBS -l mem=12GB
##PBS -o /scratch/cmz209/orthotnt/oidTest12/log/job
#PBS -o log/job
##PBS -o $OID_USER_DIR/log/job
#PBS -q cgsb-s
#PBS -N Bigmontnt
##PBS -M cmz209@nyu.edu
##PBS -m abe


echo hello for $arg1
OID_DATA=$OID_USER_DIR/data
fam=$arg1
cd $OID_DATA/$fam
echo `pwd`
echo $arg2 $arg3 $arg4
$OID_HOME/bin/runtntmx.pl "$arg2" "$arg3" "$arg4"

 

