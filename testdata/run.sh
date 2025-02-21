#!/bin/bash

export OID_USER_DIR=/scratch/$USER/testdata
export OID_HOME=~/PhyloGeneious
export ENV_WRAPPER=""

echo $OID_HOME
echo $OID_USER_DIR

nohup sh $OID_HOME/bin/topshell.sh > run.log 2>&1 &
