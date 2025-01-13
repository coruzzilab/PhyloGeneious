#!/bin/bash
#
# mod to run blastp+ by chuck zegar 10/2/2015
# slurm version 12/27/2016 Chuck Zegar
arg1=${1:-1}
arg2=${2:-1}
echo "slurm vsn called with $arg1,$arg2"
# Memory size of higher memory node needed for running mcl
HIMEM="12GB"

# Use dir name as run name
OID_RUN=$(basename $OID_USER_DIR)

JOB_SCRIPT=$OID_USER_DIR/run_test.sh
cat <<EOF >$JOB_SCRIPT
#!/bin/bash
#
# slurm job script for test blast
#


PATH=$OID_BIN:$PATH

echo "in run_test.sh $arg1"
#source activate $SOFT
#srun -n 1 -c 4  -l -o log/job/%J.out $OID_WRAPPER $OID_HOME/bin/mk_blast_parts.pl $arg1&
$OID_HOME/bin/mk_blast_parts.pl $arg1
echo hello from run_test.sh $arg1
echo slurm job \$SLURM_LOCALID
echo cpus \$SLURM_CPUS_PER_TASK
exit $JOBID

EOF
# End job script
chmod a+x $JOB_SCRIPT
#echo "$JOB_SCRIPT has $NCPU"
#            qsub -l nodes=1:ppn=$arg2 $OID_WRAPPER $JOB_SCRIPT;
#$JOB_SCRIPT
#module rm  perl
#source  deactivate
sbatch -N 1 -c 4 -o log/job/%J.out $OID_WRAPPER $JOB_SCRIPT &
sleep 2
for JOBID in $(jobs -p); do
	echo srun starts as $JOBID
done
