#!/bin/bash

args=
for i in "$@"; do 
  i="${i//\\/\\\\}"
  args="${args} \"${i//\"/\\\"}\""
done

if [ "${args}" == "" ]; then args="/bin/bash"; fi

if [[ -e /dev/nvidia0 ]]; then nv="--gpus all"; fi

#if [[ "${__HAVE_SOURCED}" == "" ]]; then
    binds="--mount type=bind,src=RUNDIR,dst=/run_folder"
#    export __HAVE_SOURCED="YES"
#fi  
    
DOCKER run ${nv} ${binds} \
phylogeneious:latest \
/bin/bash -c "
${args}
"
