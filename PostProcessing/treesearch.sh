#!/bin/bash

##72 hrs, 32 GB for mpt, jac

module purge
module load tnt/1.5
module load perl/intel/5.32.0
module load python/intel/3.8.6
#also need python 2.7

export OID_USER_DIR=/path/to/dir
export OID_HOME=/path/to/PG
export OID_BIN=$OID_HOME/bin

cd $OID_USER_DIR

# Run tree searches
echo "Tree search ..."
if [[ -f jac.tre ]]; then
	echo "Tree file already exists";
else 
	#recommend submitting mpt and jac as separate jobs
	tnt bground p $OID_HOME/PostProcessing/mpt.proc
        tnt bground p $OID_HOME/PostProcessing/jac.proc
fi

echo mpt tree done
echo jac tree done

while sleep 300; do
	if [[ ! -f mpt.tre ]]; then
		continue
	fi
	break
done

wait

#format tree results

python $OID_HOME/PostProcessing/fix_tree.py mpt.tre mpt_fixed.tre
python $OID_HOME/PostProcessing/fix_tree.py mpt.nel mpt_nel_fixed.tre
python $OID_HOME/PostProcessing/fix_tree.py jac.tre jac_fixed.tre

~/Python-2.7.16/python $OID_HOME/PostProcessing/get_node_labels.py mpt_fixed.tre mpt_nel_fixed.tre

#need mpt and jac trees
~/Python-2.7.16/python $OID_HOME/PostProcessing/parse_jac_trees.py jac_fixed.tre mpt_fixed_labeled.tre mpt_nel_labeled_with_jac_support_values.tre

python $OID_HOME/PostProcessing/labels_to_full_names.py species.txt mpt_fixed_labeled.tre mpt_nel_fixed_labeled_full_names.tre # this is the tree for PhyloBrowse/FigTree
python $OID_HOME/PostProcessing/labels_to_full_names.py species.txt mpt_nel_labeled_with_jac_support_values.tre mpt_nel_labeled_with_jac_support_values_full_names.tre # This tree will be used for generating the bpPBS.json file for PhyloBrowse (for adding the node support values)
python $OID_HOME/PostProcessing/labels_to_full_names.py species.txt patristic_distances_all.csv patristic_distances_all_full_names.csv
