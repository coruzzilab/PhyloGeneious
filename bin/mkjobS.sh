#!/bin/bash
JOB_SCRIPT=${1:-1}
cat <<EOF >$JOB_SCRIPT
#!/bin/bash
#SBATCH -o log/job/%J.out
#
# slurm job script for orthologid
#
##SBATCH -p serial
##SBATCH -n 1

export PATH=$OID_BIN:$PATH

cd \$OID_USER_DIR
MYSTART=.\$SLURM_JOB_ID.start
date +%s >\$MYSTART
echo orthologid.pl \$1 \$2 >>\$MYSTART
time
echo 'arg1 ',\$1
echo 'arg2 ',\$2

$OID_HOME/bin/orthologid.pl \$1 \$2
if [[ "\$1" == "-b" ]]; then
	#touch blast/.\$2.done
	echo ".\$2.done" >> blast/.Parts.done
fi
EOF
# End job script
