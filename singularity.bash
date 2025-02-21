#!/bin/bash

args=
for i in "$@"; do 
  i="${i//\\/\\\\}"
  args="${args} \"${i//\"/\\\"}\""
done

if [ "${args}" == "" ]; then args="/bin/bash"; fi

if [[ -e /dev/nvidia0 ]]; then nv="--nv"; fi

if [[ "${__HAVE_SOURCED}" == "" ]]; then
    binds="--bind RUNDIR:/run_folder"
    export __HAVE_SOURCED="YES"
fi  
    
SINGULARITY exec ${nv} ${binds} \
IMAGEPATH \
/bin/bash -c "
${args}
"
