#!/bin/bash
#
# mod to run blastp+ by chuck zegar 10/2/2015
# slurm version 12/27/2016 Chuck Zegar
arg1=${1:-1}
arg2=${2:-1}
arg3=${3:-1}
echo "slurm vsn called with $arg1,$arg2,$arg3"

JOB_SCRIPT=$OID_USER_DIR/run_test.sh
cat <<EOF >$JOB_SCRIPT
#!/bin/bash
#
# slurm job script for test blast
#


PATH=$OID_BIN:$PATH

echo "in run_test.sh $arg1"
#source activate $SOFT
#srun -n 1 -c 4  -l -o log/job/%J.out $OID_HOME/bin/mk_blast_parts.pl $arg1&
$OID_HOME/bin/mk_blast_partq.pl $arg1 #${ENV_WRAPPER} 
echo hello from run_test.sh $arg1
echo slurm job \$SLURM_LOCALID
echo cpus \$SLURM_CPUS_PER_TASK
exit $JOBID

EOF
# End job script
chmod a+x $JOB_SCRIPT
JOB_SCRIPTB=$OID_USER_DIR/run_blastq_job.sh
cat <<EOF >$JOB_SCRIPTB
#!/bin/bash
#SBATCH -o log/job/%J.out

echo \$SLURM_JOBID has \$SLURM_CPUS_PER_TASK cpus
$OID_HOME/bin/run_blastq.pl \$1 \$2 #${ENV_WRAPPER} 

EOF
chmod a+x $JOB_SCRIPTB
#echo "$JOB_SCRIPT has $NCPU"
#            qsub -l nodes=1:ppn=$arg2 $JOB_SCRIPT;
#$JOB_SCRIPT
#module rm  perl
#source  deactivate
if [[ ! -d $OID_USER_DIR/log/job ]]; then
	mkdir -p $OID_USER_DIR/log/job
fi
sbatch -N 1 -c $arg2 --mem $arg3 -t 04:00:00 -o log/job/%J.out $JOB_SCRIPT &
sleep 2
#for JOBID in `jobs -p`;do
#echo srun starts as $JOBID
#done
