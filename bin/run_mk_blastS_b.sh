#!/bin/bash 
#
# mod to run blastp+ by chuck zegar 10/2/2015
# slurm version 12/27/2016 Chuck Zegar
arg1=${1:-1};
arg2=${2:-1};
echo "slurm vsn called with $arg1,$arg2"
# Memory size of higher memory node needed for running mcl
HIMEM="12GB"

# Use dir name as run name
OID_RUN=$(basename $OID_USER_DIR)

JOB_SCRIPT=$OID_USER_DIR/run_test.sh
cat <<EOF >$JOB_SCRIPT
#!/bin/bash
echo "in run_test.sh $arg1"
#
# slurm job script for test blast
#

##SBATCH -p serial
#SBATCH  -o log/job/%J.out
#SBATCH  -t 120
#SBATCH   --mem $HIMEM
#SBATCH -c $arg2

PATH=$OID_BIN:$PATH

#cd \$OID_USER_DIR
#test.pl
$OID_HOME/bin/mk_blast_parts.pl $arg1

EOF
# End job script
chmod a+x $JOB_SCRIPT
echo "$JOB_SCRIPT has $NCPU"
#            qsub -l nodes=1:ppn=$arg2 $JOB_SCRIPT;
sbatch $JOB_SCRIPT&
