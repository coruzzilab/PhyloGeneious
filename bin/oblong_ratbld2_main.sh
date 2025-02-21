#!/bin/sh

#time =  00:30:00*ntask (max 08:00:00)
#to be run after all oblong_ratbld2_subprocs.sh runs have completed

family_aligned=FAMILY.aligned.revfasta

prefix=oi
fam=$1
ntimes=50
# k=0
iter=1
x=2
y=6


main() {
	out=$OID_DATADIR/$fam
	if [[ ! -f $out/oid.tre ]]; then
#merge oblong subtrees into single file
		RAND=`shuf -i1-1000000000 -n1`
		for ((i=0;i<ntimes;i++)); do 
			if [[ ! -f $out/${prefix}${i}.rat ]]; then
				echo "$out/${prefix}${i}.rat not found" 
				continue
			else
				LEN=`cat $out/${prefix}${i}.rat | wc -l`
				if [[ i -eq 0 ]]; then #first set
					head -n$((LEN-x)) $out/${prefix}${i}.rat | sed "s/);/)\*/g" > ${prefix}.rat #x y is ?
				elif [[ i -eq $((ntimes - 1)) ]]; then #last set
					tail -n$((LEN-y)) $out/${prefix}${i}.rat >> ${prefix}.rat
				else #middle sets
					head -n$((LEN-x)) $out/${prefix}${i}.rat | tail -n$((LEN-x-y)) | sed "s/);/)\*/g" >> ${prefix}.rat
				fi
			fi
		done
  		NTREE=`grep -c "^(" ${prefix}.rat`
		sed -i "s/\(tread '\)[0-9]*/\1${NTREE}/" ${prefix}.rat
		# cp ${prefix}.rat $TESTDATA/${prefix}.rat
#tree fusion with tnt
		bash $OID_HOME/bin/make_tfproc.sh $rtFamily ${prefix}.rat $out
		$ENV_WRAPPER tnt p $out/oblong_treefuse.proc
		cp oid.tre $out/oid.tre
		echo .$fam.fam.done >> $OID_DATADIR/.fam.done
#approx. tree fusion with oblong
#		NTREEMAX=$(($NTREE+10))
#		for j in {1..5}; do #damon replacememnt suggestion for tree fusion
#			# oblong -f -m5000 -i${rtFamily} -T${prefix}.rat -N -b -o${prefix}temp; #convert to nexus format
#			# cp ${prefix}temp $TESTDATA/${prefix}temp
#			cp ${prefix}.rat ${prefix}temp
#			# sed -z "s/\(\n    \)/ /g" ${prefix}temp > temp; 
 #  			if [[ $(grep -c "^(" ${prefix}temp) -gt 1 ]]; then
#				oblong -f -m5000 -i${rtFamily} -T${prefix}.rat -N -b -o${prefix}temp;
#				cp ${prefix}temp temp
#				# cp ${prefix}temp $TESTDATA/${prefix}temp
#				rm -f ${prefix}temp
#				sumtrees.py --force-rooted --min-clade-freq=0.5 temp --suppress-annotations -F newick --output-tree-filepath=${prefix}temp --no-taxa-block #get average tree
#				cp ${prefix}temp $TESTDATA/${prefix}temp
#				bash $OID_HOME/bin/nwk2tnt.sh ${prefix}temp ${prefix}tnt_temp # run nwk2tnt.sh
#			else
#				cp ${prefix}temp ${prefix}tnt_temp
#			fi
#			cp ${prefix}tnt_temp $TESTDATA/${prefix}tnt_temp
#			oblong -f -m5000 -N -i$rtFamily -T${prefix}tnt_temp -o${prefix}.tre #generate new trees
#			# cp ${prefix}.tre $TESTDATA/${prefix}.tre
#			NTREE=`grep -c "^(" ${prefix}.tre`
#			echo "oblong produced "${NTREE}" trees from sumtree"
#			if [[ $NTREE -gt $NTREEMAX ]]; then #hold
#				mv ${prefix}.tre ${prefix}temp; 
#				head -n$((7+NTREEMAX)) ${prefix}temp |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${prefix}.tre; 
#				tail -n5 ${prefix}temp >>${prefix}.tre; 
#			fi
#		done
#		cp ${prefix}.tre oids.tre
#		cp ${prefix}.tre $out/oids.tre
#		oblong -f -m5000 -i${rtFamily} -T${prefix}.tre -N -b -o$out/oid.tre;
#		NTREE=`grep -c "tree OB" $out/oid.tre`; 
#		if [[ $NTREE -gt 1 ]]; then #concatenate to single tree
#			sed -z "s/\(\n    \)/ /g" $out/oid.tre > ${prefix}temp; 
#			sumtrees.py --min-clade-freq=1.0 ${prefix}temp --suppress-annotations -F nexus --output-tree-filepath=$out/oid.tre &>> oid.log; 
#		fi
	fi
}

rtFamily=$OID_DATADIR/$fam/$family_aligned

work=$(mktemp -t -d tmp-${fam}XXX)
cd ${work}

# if [[ ${family_type[$fam]} -eq 2 ]]; then
main
# fi

rm -rf ${work}
