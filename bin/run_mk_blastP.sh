#!/bin/ksh
#
# mod to run blastp+ by chuck zegar 10/2/2015
arg1=${1:-1}
arg2=${2:-1}
print "called with $arg1,$arg2"
# Memory size of higher memory node needed for running mcl
HIMEM="12GB"

# Use dir name as run name
OID_RUN=$(basename $OID_USER_DIR)

JOB_SCRIPT=$OID_USER_DIR/run_test.sh
cat <<EOF >$JOB_SCRIPT
#!/bin/sh
#
# PBS job script for test blast
#

#PBS -S /bin/bash
#PBS -j oe
#PBS -o log/job
#PBS -l mem=$HIMEM
#PBS -l walltime=12:00:00
#PBS -N $OID_RUN
#PBS -V 

PATH=$OID_BIN:$PATH

cd \$OID_USER_DIR
#test.pl
$OID_HOME/bin/mk_blast_parts.pl $arg1 #$ENV_WRAPPER 

EOF
# End job script
chmod a+x $JOB_SCRIPT
print "$JOB_SCRIPT has $NCPU"
qsub -l nodes=1:ppn=$arg2 $JOB_SCRIPT
