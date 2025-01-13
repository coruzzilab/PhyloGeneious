export OID_USER_DIR=/scratch/$USER/testdata
export OID_HOME=~/bigplant_v4
export DATA=/scratch/$USER/testdata/PG20_data
export OID_DATADIR=/scratch/$USER/testdata/PG20_data

module load tnt/1.5
module load mcl/intel/14-137
module load perl/intel/5.32.0
module load blast+/2.13.0
module load mafft/intel/7.475
module load python/intel/3.8.6

nohup bash $OID_HOME/bin/Bigmon2.sh > $OID_USER_DIR/bigmon.log 2>&1 &