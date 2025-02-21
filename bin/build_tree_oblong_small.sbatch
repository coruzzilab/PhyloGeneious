#!/bin/bash

#SBATCH --job-name=bgm_sm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --time=12:00:00


if [ "$SLURM_ARRAY_TASK_ID" == "" ]; then exit; fi

shopt -s globstar
shopt -s dotglob nullglob

FAMILY_CNT=$OID_USER_DIR/family_cnt


num_batch=$1

source $OID_USER_DIR/family_list.sh
fams=${family_list["small"]}
read -ra families <<< "$fams"

if [[ ${ENV_WRAPPER} != "" ]]; then
    oblong_exe="$ENV_WRAPPER oblong"
else
    oblong_exe=$OID_HOME/oblong

declare -A family_type

while IFS=, read dir type || [[ $type ]]; do
  family_type["${dir}"]=${type}
done < $FAMILY_CNT

le_20 () {
	if [[ ! -f ${FAM_DIR}/oid.tre ]]; then
		ME=`pwd`
		cd ${FAM_DIR}
		$oblong_exe -f -m2000 -N -r20 -i${FILE} -ooid.tre 2>>oblong.log
		sed -n '/OPTIONS:/,$b;p' oid.tre > tempoid.tre
		$ENV_WRAPPER tnt p $OID_HOME/bin/fix.proc &> temp.log
		rm temp*
		cd ${ME}
	fi
}

le_500 () {
	if [[ ! -f ${FAM_DIR}/oid.tre ]]; then
		#drift search, 5 starting trees
		RAND=`shuf -i1-1000000000 -n1`
		timeout 7200 $oblong_exe -f -m2000 -i${FILE} -r5 -R${RAND} -d12 -o${FAM_DIR}/temp1 2>>${FAM_DIR}/oblong.log
		NTREE=`grep -c "^(" ${FAM_DIR}/temp1`; 
		if [[ ! -s ${FAM_DIR}/temp1 || $NTREE -lt 1 ]]; then
			echo "Drift error encountered: running drift search with TNT" >>${FAM_DIR}/oblong.log
			perl $OID_HOME/bin/makeDriftTree.pl ${FAM_DIR}
			return 0
		fi
		#reduces the number of 'saved' trees to improve speed - compare to TNT's 'hold 10'
		if [ $NTREE -gt 10 ]; then mv ${FAM_DIR}/temp1 temp; head -n16 temp |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${FAM_DIR}/temp1; 
		tail -n6 temp >>${FAM_DIR}/temp1; fi 

		#ratchet search
		timeout 7200 $oblong_exe -f -m2000 -i${FILE} -T${FAM_DIR}/temp1 -b -x100*4 -o${FAM_DIR}/temp2 2>>${FAM_DIR}/oblong.log
		NTREE=`grep -c "^(" ${FAM_DIR}/temp2`; 
		if [[ ! -s ${FAM_DIR}/temp2 || $NTREE -lt 1 ]]; then
			echo "Drift error encountered: running drift search with TNT" >>${FAM_DIR}/oblong.log
			# perl $OID_HOME/bin/makeDriftTree.pl ${FAM_DIR}
			return 0
		fi
		if [ $NTREE -gt 10 ]; then mv ${FAM_DIR}/temp2 temp; head -n16 temp |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${FAM_DIR}/temp2; 
		tail -n5 temp >>${FAM_DIR}/temp2; fi 

		#final search
		$oblong_exe -f -m2000 -i${FILE} -T${FAM_DIR}/temp2 -N -o${FAM_DIR}/oid.tre 2>>${FAM_DIR}/oblong.log
		

		#temporary fix for bug with current oblong version; should be able to delete later
		NTREE=`grep -c "tree O" ${FAM_DIR}/oid.tre`
		if [ $NTREE == 0 ]; then $oblong_exe -f -m2000 -i${FILE} -T${FAM_DIR}/temp2 -N -o${FAM_DIR}/oid.tre 2>>${FAM_DIR}/oblong.log; fi

		ME=`pwd`
		cd ${FAM_DIR}
		sed -n '/OPTIONS:/,$b;p' oid.tre > tempoid.tre
		$ENV_WRAPPER tnt p $OID_HOME/bin/fix.proc &> temp.log
		# rm temp*
		cd ${ME}
	fi
}

work=$(mktemp -t -d tmp-XXXXXXXX)
cd ${work}

for ((i=0;i<=num_batch;i++))
do
	index=$((SLURM_ARRAY_TASK_ID * (num_batch + 1) + i))
	echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID,num_batch=$num_batch,index=$index"
	if [[ $index -lt ${#families[*]} ]]; then
		fam=${families[$index]}
		FAM_DIR=$OID_DATADIR/$fam
		FILE=$FAM_DIR/FAMILY.aligned.revfasta
		if [[ $(grep -c ">" $FILE) -lt 4 ]];
		then
			echo "skipping fam $fam too few taxa $(grep -c ">" $FILE)";
		elif [[ ${family_type[$fam]} -eq 0 ]]; 
		then
			le_20
		elif [[ ${family_type[$fam]} -eq 1 ]];
		then
			le_500
		fi
	fi
done

rm -rf ${work}
