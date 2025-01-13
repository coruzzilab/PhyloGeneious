#!/bin/bash
#
# SLURM job script for pooltnt
#


echo hello for $1
# MYSTART=.$SLURM_JOB_ID.start
# date +%s >$MYSTART
# echo runtntmx.pl $1 $2 $3 $4 >>$MYSTART
fam=$1
#VMS 08/2022
if ! [[ -f $HOME/.passwordfile.tnt ]]; then cp $OID_HOME/.passwordfile.tnt $HOME; fi
cd $OID_DATADIR/$1
echo `pwd`
echo $2 $3 $4
perl $OID_HOME/bin/runtntmx.pl $2 $3 $4
if [[ -f $OID_DATADIR/$1/.tre.done ]]; then
	echo .$1.fam.done >> $OID_DATADIR/.fam.done
	rm $OID_DATADIR/$1/.tre.done
fi

wait
