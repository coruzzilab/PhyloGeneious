#!/bin/bash

#SBATCH --job-name=Oblong_build_tree
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2GB
#SBATCH --time=0:10:00
#SBATCH --array=0-1229
#SBATCH -o logs/job.%J.out
#SBATCH --mail-user=zt727@nyu.edu
#SBATCH --mail-type=ALL

if [ "$SLURM_ARRAY_TASK_ID" == "" ]; then exit; fi

shopt -s globstar
shopt -s dotglob nullglob

export TESTDATA=/scratch/$USER/testdata
export OBLONG_INPUT=/scratch/$USER/testdata/oblong_inputs
export OBLONG_OUTPUT=/scratch/$USER/testdata/oblong_outputs
export FAMILY_CNT=$TESTDATA/family_cnt

mkdir -p $OBLONG_OUTPUT
cd $TESTDATA

files=($OBLONG_INPUT/*)

declare -A family_type

while IFS=, read dir type || [[ $type ]]; do
  family_type["${dir}"]=${type}
done < $FAMILY_CNT

le_20 () {
	# FILE=${files[$SLURM_ARRAY_TASK_ID]} #change to input file name
	oblong -f -m2000 -N -r20 -i${FILE} -o${OBLONG_OUTPUT}/oblong_${NUMBER}.tre
}

le_500 () {
	# FILE=${files[$SLURM_ARRAY_TASK_ID]}

	#drift search, 5 starting trees
	RAND=`shuf -i1-1000000000 -n1`
	oblong -f -m2000 -i${FILE} -r5 -R${RAND} -d12 -o${OBLONG_OUTPUT}/${NUMBER}_temp1
	NTREE=`grep -c "^(" ${OBLONG_OUTPUT}/${NUMBER}_temp1`;

	temp=temp-$(cat /proc/sys/kernel/random/uuid | sed 's/[-]//g')

	#reduces the number of 'saved' trees to improve speed - compare to TNT's 'hold 10'
	if [ $NTREE -gt 10 ]; then mv ${OBLONG_OUTPUT}/${NUMBER}_temp1 ${temp}; head -n16 ${temp} |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${OBLONG_OUTPUT}/${NUMBER}_temp1; 
	tail -n6 ${temp} >>${OBLONG_OUTPUT}/${NUMBER}_temp1; fi 

	#ratchet search
	oblong -f -m2000 -i${FILE} -T${OBLONG_OUTPUT}/${NUMBER}_temp1 -b -x100*4 -o${OBLONG_OUTPUT}/${NUMBER}_temp2
	NTREE=`grep -c "^(" ${OBLONG_OUTPUT}/${NUMBER}_temp2`; 
	if [ $NTREE -gt 10 ]; then mv ${OBLONG_OUTPUT}/${NUMBER}_temp2 ${temp}; head -n16 ${temp} |sed -n "1h;2,\$H;\${g;s/\*$/;/;p}" > ${OBLONG_OUTPUT}/${NUMBER}_temp2; 
	tail -n5 ${temp} >>${OBLONG_OUTPUT}/${NUMBER}_temp2; fi 

	#final search
	oblong -f -m2000 -i${FILE} -T${OBLONG_OUTPUT}/${NUMBER}_temp2 -N -o${OBLONG_OUTPUT}/oblong_${NUMBER}.tre;

	#temporary fix for bug with current oblong version; should be able to delete later
	NTREE=`grep -c "tree O" ${OBLONG_OUTPUT}/oblong_${NUMBER}.tre`
	if [ $NTREE == 0 ]; then oblong -f -m2000 -i${FILE} -T${OBLONG_OUTPUT}/${NUMBER}_temp2 -N -o${OBLONG_OUTPUT}/oblong_${NUMBER}.tre; fi
	rm -rf ${temp}
}

# gt_500 () {
# 	#NCHAR=number of characters
# 	RAND=`shuf -i1-1000000000 -n1`
# 	#JKVAL= $NCHAR*.367 rounded up (number of jackknife positions to deactivate)

# 	#if [ -f oi${N}.tre ]; then exit; fi
# 	#if [ if oi${N} temp files ]; then _
# 	#else
# 	oblong -f -m4000 -ipal2nal.fasta -R${RAND} -ooi${N}.tmp1 
# 	#fi
# 	TMP=1
# 	#if [ $NCHAR < 3000 ]; then
# 	oblong -f -m4000 -Toi${N}.tmp1 -ipal2nal.fasta -x1 -j${JKVAL} -k -ooi${N}.tmp2 
# 	oblong -f -m4000 -Toi${N}.tmp2 -ipal2nal.fasta -ooi${N}.tmp3
# 	TMP=3
# 	#fi
# 	oblong -f -m4000 -Toi${N}.tmp${TMP} -ipal2nal.fasta -x1 -j${JKVAL} -k -ooi${N}.tmp4 
# 	oblong -f -m4000 -Toi${N}.tmp4 -ipal2nal.fasta -ooi${N}.tre
# }

if [ $SLURM_ARRAY_TASK_ID -lt ${#files[*]} ]; then
    FILE=${files[$SLURM_ARRAY_TASK_ID]}
    if [[ -f "$FILE" ]];
    then
        TMP=${FILE##*_}
        NUMBER=${TMP%.*}
        if [[ ${family_type[$NUMBER]} -eq 0 ]]; 
        then
            le_20
        elif [[ ${family_type[$NUMBER]} -eq 1 ]];
        then
            le_500
        fi
    fi
fi
