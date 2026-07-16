date
count=`ls $1/Part*faa | wc -l`
blastResDB=$2
start=1
check=0
if [[ -f $blastResDB ]]; then
	>&2 echo "Existing blastres.blst found; will continue where left off"
	GENE=`tail -n1 $blastResDB|awk '{print $1}'`
	start=`grep -w "$GENE" Part*faa |sed "s/Part\(.*\)\.faa/\1/"`
	check=1
fi
for PART in $(seq $start $count); do
	[[ ! -f $1/Part${PART}.faa ]] && >&2 echo "Part file Part${PART}.faa not found!";
	for GENE in `grep ">" $1/Part${PART}.faa`; do
		if [[ $check > 0 ]]; then
			if [[ $(grep -m1 -w "$GENE" $blastResDB|wc -l) > 0 ]]; then continue; else check=0; fi
		fi
		echo $GENE >> $blastResDB
		grep -w ^${GENE:1} */*Part${PART}.blst |awk '!a[$1$2]++'|cut -d":" -f2- >> $blastResDB
	done
	check=0
done 
wait
date
