#!/bin/bash
#
# SLURM job script for pooltnt
#
#SBATCH -o log/job/tntmx.%J.out

echo hello for $1
OID_DATA=$OID_USER_DIR/data
MYSTART=.$SLURM_JOB_ID.start
date +%s >$MYSTART
echo runtntmx.pl $1 $2 $3 $4 >>$MYSTART
fam=$1
#VMS 08/2022
if ! [[ -f $HOME/.passwordfile.tnt ]]; then cp $OID_HOME/.passwordfile.tnt $HOME; fi
cd $OID_DATA/$fam
echo $(pwd)
echo $2 $3 $4
$OID_HOME/bin/runtntmx.pl $2 $3 $4
