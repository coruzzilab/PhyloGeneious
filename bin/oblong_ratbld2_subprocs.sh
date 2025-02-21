
family_aligned=FAMILY.aligned.revfasta

prefix=oi
fam=$1
kk=$2
# k=0
iter=1

if [[ ${ENV_WRAPPER} != "" ]]; then
    oblong_exe="$ENV_WRAPPER oblong"
else
    oblong_exe=$OID_HOME/oblong

# declare -A family_type

# while IFS=, read dir type || [[ $type ]]; do
#   family_type["${dir}"]=${type}
# done < $FAMILY_CNT

le_1500 () {
	out=$OID_DATADIR/$fam
	if [[ ! -f $out/oid.tre ]]; then
		for ((j=0;j<5;j++)); do
			k=$((kk * 5 + j))
			if [[ ! -f $out/${prefix}${k}.rat ]]; then
				RAND=`shuf -i1-1000000000 -n1`
				#collapse rule 3?
				#. as character?
				$oblong_exe -f -m4000 -R${RAND} -i${rtFamily} -o${prefix}${k}temp1 2>>${out}/oblong.log;
				NTREE=`grep -c "^(" ${prefix}${k}temp1`; 
				if [[ $NTREE -gt 1 ]]; then #hold 1
					mv ${prefix}${k}temp1 ${prefix}${k}temp; 
					head -n7 ${prefix}${k}temp |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${prefix}${k}temp1; 
					tail -n6 ${prefix}${k}temp >>${prefix}${k}temp1; 
				fi
				for ((i=1;i<=iter;i++)); do
					if [[ -f ${prefix}${k}temp3 ]]; then
						$oblong_exe -f -m4000 -i${rtFamily} -T${prefix}${k}temp3 -x1*10 -j37 -o${prefix}${k}temp2 2>>${out}/oblong.log;
					else
						$oblong_exe -f -m4000 -i${rtFamily} -T${prefix}${k}temp1 -x1*10 -j37 -o${prefix}${k}temp2 2>>${out}/oblong.log;
					fi
					NTREE=`grep -c "^(" ${prefix}${k}temp2`; 
					if [[ $NTREE -gt 2 ]]; then #hold 2
						mv ${prefix}${k}temp2 ${prefix}${k}temp; 
						head -n8 ${prefix}${k}temp |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${prefix}${k}temp2; 
						tail -n5 ${prefix}${k}temp >>${prefix}${k}temp2; 
					fi
					$oblong_exe -f -m4000 -i${rtFamily} -T${prefix}${k}temp2 -o${prefix}${k}temp3 2>>${out}/oblong.log;
					NTREE=`grep -c "^(" ${prefix}${k}temp3`; 
					if [[ $NTREE -gt 2 && i -lt $iter ]]; then #hold 2
						mv ${prefix}${k}temp3 ${prefix}${k}temp; 
						head -n8 ${prefix}${k}temp |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${prefix}${k}temp3; 
						tail -n5 ${prefix}${k}temp >>${prefix}${k}temp3; 
					fi
				done
				#ko?
				if [[ $NTREE -gt $(($iter*3 + 1)) ]]; then #hold 
					mv ${prefix}${k}temp3 ${prefix}${k}temp; 
					head -n${7 + $iter*3} ${prefix}${k}temp |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${prefix}${k}temp3; 
					tail -n5 ${prefix}${k}temp >>${prefix}${k}temp3; 
				fi
				sed -n '/OPTIONS:/,$b;p' ${prefix}${k}temp3 > $out/${prefix}${k}.rat
				# oblong -f -m4000 -i${rtFamily} -T${prefix}${k}temp3 -N -b -o$out/${prefix}${k}.tre;
			fi
			# NTREE=`grep -c "tree OB" $out/${prefix}${k}.tre`; 
			# if [[ $NTREE -gt 1 ]]; then #concatenate to single tree
			# 	sed -z "s/\(\n    \)/ /g" $out/${prefix}${k}.tre > ${prefix}${k}temp; 
			# 	sumtrees.py --min-clade-freq=1.0 ${prefix}${k}temp --suppress-annotations -F nexus --output-tree-filepath=$out/${prefix}${k}.tre &>> ${prefix}${k}.log; 
			# fi
		done
	fi

}


# fam=47
rtFamily=$OID_DATADIR/$fam/$family_aligned

work=$(mktemp -t -d tmp-${fam}XXXX)
cd ${work}

# if [[ ${family_type[$fam]} -eq 2 ]]; then
le_1500
# fi

rm -rf ${work}