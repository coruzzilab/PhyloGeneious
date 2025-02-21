#!/bin/bash
#
# PBS job script for orthologid
#

#PBS -V 
#PBS -S /bin/bash
#PBS -j oe
#PBS -o log/job
#PBS -N Bigmontnt

echo hello for $arg1
OID_DATA=$OID_USER_DIR/data
fam=$arg1
if ! [[ -f $HOME/.passwordfile.tnt ]]; then cp $OID_HOME/.passwordfile.tnt $HOME; fi
cd $OID_DATA/$fam
echo `pwd`
echo $arg2 $arg3 $arg4
$OID_HOME/bin/runtntmx.pl "$arg2" "$arg3" "$arg4"  #$ENV_WRAPPER 


