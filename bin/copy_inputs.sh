#!/bin/bash

export TRANSLATED_DIR="/scratch/zt727/testdata/PG20_reverse_translated"
export DATA_DIR="/scratch/zt727/testdata/data_parallel"

cnt=0
for dir in "$TRANSLATED_DIR"/*
do
    num=${dir##*/}
    mv $TRANSLATED_DIR/$num/FAMILY.aligned $TRANSLATED_DIR/$num/FAMILY.aligned.revfasta
    cp $DATA_DIR/$num/oid.nex $TRANSLATED_DIR/$num/
    cp $DATA_DIR/$num/FAMILY $TRANSLATED_DIR/$num/
    cp $DATA_DIR/$num/FAMILY.aligned $TRANSLATED_DIR/$num/
    let cnt++
done

echo "$cnt families copied"