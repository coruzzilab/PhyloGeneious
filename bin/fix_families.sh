export DATA=/scratch/$USER/testdata/PG20_data
export OID_USER_DIR=/scratch/$USER/testdata
export OID_HOME=~/bigplant_v4


readarray -t families < $OID_USER_DIR/FAMILY

for fam in "${families[@]}"
do
    bash $OID_HOME/bin/fix_familyfastas.sh $DATA/$fam/FAMILY
done