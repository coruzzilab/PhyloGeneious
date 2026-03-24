#!/bin/bash

module load tnt/1.5
module load python/intel/3.8.6

export OID_USER_DIR=/path/to/dir
export OID_HOME=/path/to/PG
export ASTRAL=/path/to/astral-pro3

#get input gene trees file with branchlengths
for DIR in `ls $OID_USER_DIR/data/*/oid.tre | sed "s/oid.tre//"`; do 
cd $DIR; 
if [[ ! -s brlen.tre ]]; then 
	sed "s/#/_/g" oid.nex > brtemp.nex; 
	sed "s/#/_/g" oid.tre > brtemp.tre; 
	sed -i "s/_NEXUS/#NEXUS/" brtemp.tre; 
	tnt p $OID_HOME/PostProcessing/getbrlengths.proc; 
	python $OID_HOME/PostProcessing/fix_tree.py brtemp2.tre brlen.tre; 
	sed -i "s/_[^,();: ]*//g" brlen.tre 
	cat brlen.tre >> $OID_USER_DIR/all_mp_famtrees_br.nw; 
	rm brtemp*
fi
cd $OID_USER_DIR; 
done

#run ASTRALPro
$ASTRAL -i $OID_USER_DIR/all_mp_famtrees_br.nw -o $OID_USER_DIR/astralpro.tre 2>$OID_USER_DIR/astralpro.log
