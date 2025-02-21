prefix=pi
family_aligned=FAMILY.aligned.revfasta
LOGS=$OID_USER_DIR/log

source $OID_USER_DIR/family_list.sh
fams=${family_list["ratbld6"]}
read -ra families <<< "$fams"


main() {
if [[ ! -f $OID_DATADIR/$fam/oid.tre ]]; then
  bash $OID_HOME/bin/setup_proc_dna.sh $fam
# if [[ "$HPC" == "S" ]] ; then
  sbatch --job-name=bgm_alt6 -n 1 -c 4 -t 96:00:00 --mem 128GB -o $LOGS/oblong_ratbld6_alternate_logs/pi.%J.out $OID_HOME/bin/run_tntmx_dna${HPC}.sh $fam 8 4 pi;
fi
# bash $OID_HOME/bin/run_tntmx_dna${HPC}.sh $fam 8 4 pi;
#run normal TNT programming
#arguments: 8=number subtrees to make; 4=number of threads; pi=prefix
# else
# qsub -l nodes=1:ppn=$NCPU -l mem=32GB,walltime=96:00:00 $OID_WRAPPER run_tntmx_dna${HPC}.sh -v arg1 $fam arg2 8 arg3 4 arg4 pi;
# fi
}


for fam in "${families[@]}"
do
  #fam=${families[$i]}
  # if [[ ${family_type[$fam]} -eq 6 ]]; then
  main
  # fi
done

