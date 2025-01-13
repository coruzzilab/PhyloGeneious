#!/bin/ksh
MY_SHELL=$0
export MY_SHELL
OID_USER_DIR=$(pwd)
export OID_USER_DIR
echo $OID_HOME
echo $OID_USER_DIR
cd $OID_USER_DIR
qsub $OID_HOME/bin/pipe.pbs
