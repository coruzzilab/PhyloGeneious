#!/bin/bash

export OID_USER_DIR=/scratch/$USER/testdata
export OID_HOME=/~
export ENV_WRAPPER=""

#module load blast+/x86_64/2.11.0+
#module load mcl/x86_64/22-282
#module load mafft/x86_64/7.453
#module load diamond/x86_64/2.1.11
#module load MMseqs2/x86_64/14-7e284

echo $OID_HOME
echo $OID_USER_DIR

nohup sh $OID_HOME/bin/topshell.sh > run.log 2>&1 
tail -n1 run.log
