#!/bin/bash
#
# PBS job script for orthologid
#

#PBS -V 
#PBS -S /bin/bash
#PBS -j oe
#PBS -o log/job
##PBS -o $OID_USER_DIR/log/job
#PBS -q cgsb-s
#PBS -N Bigmontnt

echo hello for $arg1
OID_DATA=$OID_USER_DIR/data
cd $OID_USER_DIR
MYSTART=.$PBS_JOBID.start
$OID_HOME/bin/gettime.pl >$MYSTART
echo runtntmx.pl  "$arg1" "$arg2" "$arg3" "$arg4" >>$MYSTART
if ! [[ -f $HOME/.passwordfile.tnt ]]; then cp $OID_HOME/.passwordfile.tnt $HOME; fi
fam=$arg1
cd $OID_DATA/$fam
echo `pwd`
echo $arg2 $arg3 $arg4
$OID_HOME/bin/runtntmx.pl "$arg2" "$arg3" "$arg4"

 

