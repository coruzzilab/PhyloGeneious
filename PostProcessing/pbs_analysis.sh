#!/bin/bash

module load tnt/1.5
module load perl/intel/5.32.0

export OID_USER_DIR=/path/to/dir
export OID_HOME=/path/to/PG
export HPC=S

### run PBS
mkdir $OID_USER_DIR/pbscompare
cd $OID_USER_DIR/pbscompare
perl $OID_HOME/PostProcessing/pbs_split_bowery_gil.pl -l pbs.log -m $OID_USER_DIR/Matrix.tnt -n $OID_USER_DIR/mpt.tre -t $OID_USER_DIR/mpt.tre > pbs.csv

