#!/bin/sh
export DATA=$OID_USER_DIR/data


readarray -t families < $OID_USER_DIR/FAMILY

for fam in "${families[@]}"
do
    bash $OID_HOME/bin/fix_familyfastas.sh $DATA/$fam/FAMILY
done
