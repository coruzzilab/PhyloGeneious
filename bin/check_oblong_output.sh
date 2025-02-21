OID_DATADIR=$OID_USER_DIR/data
source $OID_USER_DIR/family_list.sh
echo "checking oblong output trees"
missing=false
for type in "${!family_list[@]}"; do
    echo "=== checking output trees for type $type families ==="
    fams=${family_list[$type]}
    read -ra fam_list <<< "$fams"
    missing_output=()
    for fam in "${fam_list[@]}"
    do
        if [[ ! -s $OID_DATADIR/$fam/oid.tre ]]; then
            missing_output+=($fam)
        else
            NTREE=`grep -c "^(" $OID_DATADIR/$fam/oid.tre`;
            if [[ $NTREE -lt 1 ]]; then
                missing_output+=($fam)
            fi
        fi
    done
    if [[ ${#missing_output[@]} > 0 ]]; then
        echo "type $type output tree missing, ${#missing_output[@]} missing families: ${missing_output[@]}"
        missing=true
    else
        echo "type $type done"
    fi
done