#!/bin/bash
#
# PBS job script for orthologid
#

#PBS -S /bin/bash
#PBS -j oe
#PBS -o log/job
#PBS -l mem=12GB
##PBS -o /scratch/cmz209/orthotnt/oidTest9/log/job/
#PBS -q cgsb-s
#PBS -N pooltnt
#PBS -V

cd $OID_USER_DIR
MYSTART=.$PBS_JOBID.start
date +%s >$MYSTART
echo pooltnt.pl "$arg1" "$arg2" "$arg3" >>$MYSTART
date
time
$OID_HOME/bin/pooltnt.pl "$arg1" "$arg2" "$arg3"

date
time
