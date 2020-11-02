#!/bin/bash
#
# PBS job script for orthologid
#

#PBS -S /bin/bash
#PBS -j oe
#PBS -l mem=$HIMEM
##PBS -o /scratch/cmz209/orthotnt/oidTest9/log/job/
#PBS -o log/job/
#PBS -q 
#PBS -N 
#PBS -V 

PATH=:/usr/local/bin:/Users/chuckza/bin:/Users/chuckza/perlbin:/Users/chuckza/perl5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/usr/local/git/bin:/usr/local/ncbi/blast/bin:/usr/texbin:/Users/chuckza/blastplus/bin

cd $OID_USER_DIR
MYSTART=.$PBS_JOBID.start
/Users/chuckza/bigplant_v4/bin/gettime.pl >$MYSTART
echo orthologid.pl $1 $2 >>$MYSTART
date
time
orthologid.pl "$arg1" "$arg2"

date
time
if [[ "$arg1" == "-b" ]]; then
	touch blast/.$arg2.done
fi
