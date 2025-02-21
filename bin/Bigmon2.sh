#!/bin/bash

# set -e

export OID_DATADIR=$OID_USER_DIR/data
LOGS=$OID_USER_DIR/log

# get families
families_dir=($OID_DATADIR/*)
raw_families=("${families_dir[@]##*/}")

families=()
for family in "${raw_families[@]}"
do
    case $family in ''|*[!0-9]*) 
        continue;;
    esac;
    families+=($family)
done

echo "num of families: ${#families[@]}"
printf "%s\n" "${families[@]}" > $OID_USER_DIR/FAMILY

# job config
subprocs_ntimes=10
max_jobs=$1
num_batch=1
success=1

arrayqs () {
    script=$1
    outlog=$2
    lastarray=$(($3 - 1))
    args=$4
    
    if [[ $HPC -eq "S" ]]; then
        sbatch --array=0-$lastarray -o $LOGS/$outlog $script $args
    else
        qsub -t 0-$lastarray -o $LOGS/mafft_logs/%J.out $script $args
    fi
}

getjobs () {
    if [[ $HPC -eq "S" ]]; then
        squeue -u $USER | tail -n +2 
    else
        qstat -e -u $MYUSER #check USER variable
    fi
    success=$?
}

# mafft

echo "$(date +"%Y-%m-%d %H:%M:%S") start running mafft"
num_jobs=${#families[@]}
num_batch=$(($num_jobs / $max_jobs))
if [[ $num_batch > 0 ]]; then
    num_jobs=$max_jobs
fi
rm -rf $LOGS/mafft_logs
mkdir -p $LOGS/mafft_logs
arrayqs $OID_HOME/bin/run_mafftS.sbatch mafft_logs/%J.out $num_jobs $num_batch
#sbatch --array=0-$(($num_jobs - 1)) -o $LOGS/mafft_logs/%J.out $OID_HOME/bin/run_mafftS.sbatch $num_batch

running_jobs=$(getjobs | grep "run_maff")

while [[ ! -z "$running_jobs" ]];
do
    sleep 1m
    running_jobs=$(getjobs | grep "run_maff")
done
echo "$(date +"%Y-%m-%d %H:%M:%S") mafft jobs finished"

echo "checking mafft output files"
missing_output=()
for fam in "${families[@]}"
do
    if [[ ! -f $OID_DATADIR/$fam/oid.nex ]]; then
        missing_output+=($fam)
    fi
done
if [[ ${#missing_output[@]} > 0 ]]; then
    echo "type $type mafft output missing, ${#missing_output[@]} missing families: ${missing_output[@]}"
    echo "aborting"
    exit 1
fi
echo "mafft completed"

# count family size
echo "counting family size"
$ENV_WRAPPER python $OID_HOME/bin/count_families.py $OID_DATADIR $OID_USER_DIR/family_cnt

declare -A family_type
while IFS=, read dir type || [[ $type ]]; do
  family_type["${dir}"]=${type}
done < $OID_USER_DIR/family_cnt

# invert the dict from {family -> type} to {type -> list of families}
declare -A type_map=( [0]="small" [1]="small" [2]="ratbld2" [3]="ratbld3" [4]="ratbld3" [5]="ratbld3" [6]="ratbld6" [7]="toobig")
declare -A family_list
for fam in "${!family_type[@]}"; do
    type=${type_map[${family_type[$fam]}]}
    if [ -z "${family_list[$type]}" ]; then
        family_list[$type]="$fam"
    else
        family_list[$type]="${family_list[$type]} $fam"
    fi
done
declare -p family_list > $OID_USER_DIR/family_list.sh

# declare -A family_cnt

# reverse translation
echo "start reverse translating"
$ENV_WRAPPER python $OID_HOME/bin/reverse_translate_families.py $OID_DATADIR $OID_DATADIR

echo "checking reverse translate output files"
missing_output=()
for fam in "${families[@]}"
do
    if [[ ! -f $OID_DATADIR/$fam/FAMILY.aligned.revfasta ]]; then
        missing_output+=($fam)
    fi
done
if [[ ${#missing_output[@]} > 0 ]]; then
    echo "reverse translate output files missing, missing families: ${missing_output[@]}"
    echo "aborting"
    exit 1
fi
echo "reverse translate completed"

# oblong tree generation

# track job status
declare -A job_status=( ["small"]="INIT" ["ratbld2"]="INIT" ["ratbld3"]="INIT" ["ratbld6"]="INIT")

for type in "${!family_list[@]}"; do
    echo "=== generate trees for type $type families ==="
    fams="${family_list[$type]}"
    read -ra fam_list <<< "$fams"
    num_fams=${#fam_list[@]}
    if [[ $type = "small" ]]; then
        prefix=oi
        num_jobs=$num_fams
        # echo "num_jobs=$num_jobs"
        num_batch=$(($num_jobs / $max_jobs))
        if [[ $num_batch > 0 ]]; then
            num_jobs=$max_jobs
        fi
        echo "num_batch=$num_batch"
        rm -rf $LOGS/oblong_small_logs
        mkdir -p $LOGS/oblong_small_logs
        arrayqs $OID_HOME/bin/build_tree_oblong_small.sbatch oblong_small_logs/${prefix}.%J.out $num_jobs $num_batch
        #sbatch --array=0-$(($num_jobs - 1)) -o $LOGS/oblong_small_logs/${prefix}.%J.out $OID_HOME/bin/build_tree_oblong_small.sbatch $DATA $num_batch
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [0,1]: INIT -> RUNNING"
        job_status["small"]="RUNNING"
    elif [[ $type = "ratbld2" ]]; then
        prefix=oi
        num_jobs=$(($subprocs_ntimes * $num_fams))
        num_batch=$(($num_jobs / $max_jobs))
        if [[ $num_batch > 0 ]]; then
            num_jobs=$max_jobs
        fi
        rm -rf $LOGS/oblong_ratbld2_subproc_logs
        mkdir -p $LOGS/oblong_ratbld2_subproc_logs
        arrayqs $OID_HOME/bin/oblong_ratbld2_subprocs.sbatch oblong_ratbld2_subproc_logs/${prefix}.%J.out $num_jobs $num_batch
        #sbatch --array=0-$(($num_jobs - 1)) -o $LOGS/oblong_ratbld2_subproc_logs/${prefix}.%J.out $OID_HOME/bin/oblong_ratbld2_subprocs.sbatch $DATA $num_batch
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [2]: INIT -> SUBPROCS_RUNNING"
        job_status["ratbld2"]="SUBPROCS_RUNNING"
    elif [[ $type = "ratbld3" ]]; then
        prefix=op
        num_jobs=$(($subprocs_ntimes * $num_fams))
        num_batch=$(($num_jobs / $max_jobs))
        if [[ $num_batch > 0 ]]; then
            num_jobs=$max_jobs
        fi
        rm -rf $LOGS/oblong_ratbld3_subproc_logs
        mkdir -p $LOGS/oblong_ratbld3_subproc_logs
        arrayqs $OID_HOME/bin/oblong_ratbld3_subprocs.sbatch oblong_ratbld3_subproc_logs/${prefix}.%J.out $num_jobs $num_batch
        #sbatch --array=0-$(($num_jobs - 1)) -o $LOGS/oblong_ratbld3_subproc_logs/${prefix}.%J.out $OID_HOME/bin/oblong_ratbld3_subprocs.sbatch $DATA $num_batch
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [3,5]: INIT -> SUBPROCS_RUNNING"
        job_status["ratbld3"]="SUBPROCS_RUNNING"
    elif [[ $type = "ratbld6" ]]; then
        prefix=pi
        rm -rf $LOGS/oblong_ratbld6_alternate_logs
        mkdir -p $LOGS/oblong_ratbld6_alternate_logs
        bash $OID_HOME/bin/oblong_ratbld6_alternate.sh
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [6]: INIT -> RUNNING"
        job_status["ratbld6"]="RUNNING"
    fi
done

#run small jobs in serial - save job resources
#echo "=== generate trees for very small families ==="
# for type in "${!family_list[@]}"; do
#     fams="${family_list[$type]}"
#     read -ra fam_list <<< "$fams"
#     num_fams=${#fam_list[@]}
#     if [[ $type = "vsmall" ]]; then #run serial - saves job resources
#         for fam in ${fam_list[@]}
#         do
#             FAM_DIR=$OID_DATADIR/$fam
#             FILE=$FAM_DIR/FAMILY.aligned.revfasta
#             if [[ $(grep -c ">" $FILE) -lt 4 ]]; then
#                 echo "skipping fam $fam too few taxa $(grep -c ">" $FILE)";
#             else
#                 if [[ ! -f ${FAM_DIR}/oid.tre ]]; then
#                     oblong -f -m2000 -N -r20 -i${FILE} -o${FAM_DIR}/oid.tre 2>>${FAM_DIR}/oblong.log
#                     ME=`pwd`
#                     cd ${FAM_DIR}
#                     sed -n '/OPTIONS:/,$b;p' oid.tre > tempoid.tre
#                     tnt p $OID_HOME/bin/fix.proc &> temp.log
#                     rm temp*
#                     cd ${ME}
#                 fi
#             fi
#         done
#     fi
# done

# check job status and advance

while true;
do
    mapfile -t job_list < <( getjobs | grep "bgm" ) 
    # if [[ $success == 0 ]]; then
    unset active_job
    declare -A active_job
    for job in "${job_list[@]}"
    do
        curr_job=$(echo $job | awk '{print $3}')
        active_job[$curr_job]=1
    done

    if [[ ${job_status["ratbld2"]} = "MAIN_RUNNING" ]] && [[ ! ${active_job["bgm_m2"]} ]]; then
        echo "ratbld2 main finished"
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [2]: MAIN_RUNNING -> DONE"
        job_status["ratbld2"]="DONE"
    fi
    if [[ ${job_status["ratbld3"]} = "MAIN_RUNNING" ]] && [[ ! ${active_job["bgm_m3"]} ]]; then
        echo "ratbld3 main finished"
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [3,5]: MAIN_RUNNING -> DONE"
        job_status["ratbld3"]="DONE"
    fi
    if [[ ${job_status["small"]} = "RUNNING" ]] && [[ ! ${active_job["bgm_sm"]} ]]; then
        echo "small finished"
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [0,1]: RUNNING -> DONE"
        job_status["small"]="DONE"
    fi
    if [[ ${job_status["ratbld6"]} = "RUNNING" ]] && [[ ! ${active_job["bgm_alt6"]} ]]; then
        echo "ratbld6 finished"
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [6]: RUNNING -> DONE"
        job_status["ratbld6"]="DONE"
    fi

    if [[ ${job_status["ratbld2"]} = "SUBPROCS_RUNNING" ]] && [[ ! ${active_job["bgm_s2"]} ]]; then
        echo "=== ratbld2 subprocs finished, starting ratbld2 main ==="
        prefix=oi
        fams=${family_list["ratbld2"]}
        read -ra fam_list <<< "$fams"
        num_fams=${#fam_list[@]}
        num_jobs=$num_fams
        num_batch=$(($num_jobs / $max_jobs))
        if [[ $num_batch > 0 ]]; then
            num_jobs=$max_jobs
        fi
        rm -rf $LOGS/oblong_ratbld2_main_logs
        mkdir -p $LOGS/oblong_ratbld2_main_logs
        arrayqs $OID_HOME/bin/oblong_ratbld2_main.sbatch oblong_ratbld2_main_logs/${prefix}.%J.out $num_jobs $num_batch
        #sbatch --array=0-$(($num_jobs - 1)) -o $LOGS/oblong_ratbld2_main_logs/${prefix}.%J.out $OID_HOME/bin/oblong_ratbld2_main.sbatch $DATA $num_batch
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [2]: SUBPROCS_RUNNING -> MAIN_RUNNING"
        job_status["ratbld2"]="MAIN_RUNNING"
    fi
    if [[ ${job_status["ratbld3"]} = "SUBPROCS_RUNNING" ]] && [[ ! ${active_job["bgm_s3"]} ]]; then
        echo "=== ratbld3 subprocs finished, starting ratbld3 main ==="
        prefix=op
        fams=${family_list["ratbld3"]}
        read -ra fam_list <<< "$fams"
        num_fams=${#fam_list[@]}
        num_jobs=$num_fams
        num_batch=$(($num_jobs / $max_jobs))
        if [[ $num_batch > 0 ]]; then
            num_jobs=$max_jobs
        fi
        rm -rf $LOGS/oblong_ratbld3_main_logs
        mkdir -p $LOGS/oblong_ratbld3_main_logs
        arrayqs $OID_HOME/bin/oblong_ratbld3_main.sbatch oblong_ratbld3_main_logs/${prefix}.%J.out $num_jobs $num_batch
        #sbatch --array=0-$(($num_jobs - 1)) -o $LOGS/oblong_ratbld3_main_logs/${prefix}.%J.out $OID_HOME/bin/oblong_ratbld3_main.sbatch $DATA $num_batch
        echo "$(date +"%Y-%m-%d %H:%M:%S") [oblong] type [3,5]: SUBPROCS_RUNNING -> MAIN_RUNNING"
        job_status["ratbld3"]="MAIN_RUNNING"
    fi

    declare -p job_status > $LOGS/job_status
    mapfile -t job_list < <( getjobs | grep "bgm" )
    if [[ ${#job_list[@]} = 0 ]]; then
        break
    fi
    sleep 5m
    # else
    #     echo "queue query failed"
    #     sleep 30
    # fi
done
# check oblong output files
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

echo "oblong tree generation completed"
echo "$(date +"%Y-%m-%d %H:%M:%S") Bigmon Done!"
>$LOGS/job/bigmdone
