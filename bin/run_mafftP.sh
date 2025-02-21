#!/bin/bash
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
MYSTART=.$PBS_JOBID.start
date +%s >$MYSTART
echo mafft "$arg1" "$arg2" >>$MYSTART

cd $OID_DATA/"$arg1"
$ENV_WRAPPER mafft --auto --quiet --anysymbol --thread "$arg2" FAMILY >FAMILY.aligned
if [[ ! -s FAMILY.aligned ]]; then
	echo "failed re-run mafft "
	$ENV_WRAPPER mafft --nofft --retree 1 --memsavetree --quiet --anysymbol --thread $arg2 FAMILY >FAMILY.aligned
fi
cd $OID_USER_DIR
$OID_HOME/bin/makenex.pl "$arg1" #$ENV_WRAPPER 

date
time
