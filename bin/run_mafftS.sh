#!/bin/bash
#
# SLURM job script for pooltnt
#
#SBATCH -o log/job/mafft.%J.out


cd $OID_USER_DIR
date
time

echo hello for mafft
OID_DATA=$OID_USER_DIR/data
MYSTART=.$SLURM_JOB_ID.start
$OID_HOME/bin/gettime.pl >$MYSTART
echo mafft $1 $2  >>$MYSTART

cd $OID_DATA/$1
mafft --quiet --anysymbol --thread $2 FAMILY>FAMILY.aligned
#wc FAMILY.aligned
if [[ ! -s FAMILY.aligned ]]; then
    echo "failed re-run mafft "
    mafft --nofft --retree 1 --memsavetree --quiet --anysymbol --thread $2 FAMILY>FAMILY.aligned
fi
cd $OID_USER_DIR
$OID_HOME/bin/makenex.pl $1 


date
time
