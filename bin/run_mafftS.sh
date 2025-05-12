#!/bin/bash
#
#SBATCH -o log/job/mafft.%J.out

cd $OID_USER_DIR
date
time

echo hello for mafft
if [[ $OID_DATADIR ]]; then 
	OID_DATA=$OID_DATADIR
else OID_DATA=$OID_USER_DIR/data
fi
# MYSTART=.$SLURM_JOB_ID.start
# date +%s >$MYSTART
# echo mafft $1 $2 >>$MYSTART

cd $OID_DATA/$1

if [[ ! -f oid.nex ]]; then
	$ENV_WRAPPER mafft --quiet --anysymbol --thread $2 FAMILY >FAMILY.aligned
	#wc FAMILY.aligned
	if [[ ! -s FAMILY.aligned ]]; then
		echo "failed re-run mafft "
		$ENV_WRAPPER mafft --nofft --retree 1 --memsavetree --quiet --anysymbol --thread $2 FAMILY >FAMILY.aligned
	fi
	cd $OID_USER_DIR
	$OID_HOME/bin/makenex.pl $1 #$ENV_WRAPPER 
fi
date
time
