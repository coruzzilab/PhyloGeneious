#!/bin/bash
JOB_SCRIPT=${1:-1}
# Generate job script
cat <<EOF >$JOB_SCRIPT
#!/bin/bash
#
# PBS job script for orthologid
#

#PBS -S /bin/bash
#PBS -j oe
#PBS -l mem=\$HIMEM
#PBS -o log/job/
#PBS -N $OID_RUN
#PBS -V

#PATH=$OID_BIN:$PATH

cd \$OID_USER_DIR
MYSTART=.\$PBS_JOBID.start
date +%s >\$MYSTART
echo orthologid.pl "\$arg1" "\$arg2" >>\$MYSTART
date
time
$OID_HOME/bin/orthologid.pl "\$arg1" "\$arg2" #$ENV_WRAPPER 

date
time
if [[ "\$arg1" == "-b" ]]; then
	#touch blast/.\$arg2.done
	echo ".\$arg2.done" >> blast/.Parts.done
fi
EOF
# End job script
