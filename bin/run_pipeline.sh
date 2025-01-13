#!/bin/bash
#
#
# Script for running OrthologID pipeline
# This script makes use of PBS for job submission.
#
# Copyright (C) 2006-2012 Ernest K. Lee
#
#    This file is part of OrthologID.
#
#    OrthologID is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    OrthologID is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with OrthologID.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Author: Ernest K Lee <elee@amnh.org>
#
# mod to run blastp+ by chuck zegar 10/2/2015
# My PBS queue
PBSQ="cgsb-s"

# Memory size of higher memory node needed for running mcl
HIMEM="128GB"
##export OID_USER_DIR=`pwd`
##export OID_HOME=/home/cmz209/orthotnt/OID_nw3
if [[ -z $MY_SHELL ]]; then
	echo "MY_SHELL not defined - won't restart properly"
	export MY_SHELL=$0
##	exit 1
fi
echo "will restart $MY_SHELL"
date +%s >starttime
#module load mercurial/intel/2.1.2
#module load mcl/intel/12-068
#module load mafft/intel/6.864

# Check environment variables and user directory
#if ! env | grep -q "^OID_HOME="; then
#	echo  "OID_HOME not defined ... exiting"
#	exit 1
#fi

# Get version
if [[ -f $OID_HOME/VERSION ]]; then
	OID_VERSION=$(<$OID_HOME/VERSION)
else
	OID_VERSION="unknown"
fi
if [[ ! -f $OID_USER_DIR/procfiles.txt ]]; then
	cp $OID_HOME/config/procfiles.txt $OID_USER_DIR
	echo "copied procfiles.txt from $OID_HOME"
fi

OID_BIN=$OID_HOME/bin
echo "oidbin $OID_BIN"

#if [[ $# -eq 0 ]] && ! env | grep -q "^OID_USER_DIR="; then
#	echo  "Must specify run directory as argument or define OID_USER_DIR"
#	exit 1
#elif [[ $# -eq 1 ]]; then
#	export OID_USER_DIR="$1"
#elif [[ $# -gt 1 ]]; then
#	echo  "Usage: $0 [ OID_USER_DIR ]"
#	exit 2
#fi
#if [[ ! -d $OID_USER_DIR ]]; then
#	echo  "OID_USER_DIR \"$OID_USER_DIR\" is not a directory"
#	exit 1
#fi

# Use dir name as run name
OID_RUN=$(basename $OID_USER_DIR)

echo "== Starting OrthologID version $OID_VERSION pipeline =="
echo "Run directory is $OID_USER_DIR"

# Check/parse config file
CONFIG=$OID_USER_DIR/config
if [[ ! -f $CONFIG ]]; then
	echo "Config file required -- a sample is $OID_HOME/config/config "
	echo "... exiting"
	exit 1
fi
declare -a INGROUP
INGROUP=($(grep INGROUP $CONFIG | cut -f2 -d=))
declare -a OUTGROUP
OUTGROUP=($(grep OUTGROUP $CONFIG | cut -f2 -d=))
#set -a INGROUP $(grep INGROUP $CONFIG | cut -f2 -d=)
#set -a OUTGROUP $(grep OUTGROUP $CONFIG | cut -f2 -d=)
HPC=$(sed -n 's/^HPC *= *\([PSG]*\).*/\1/p' $CONFIG)
if [[ -z $HPC ]]; then
	HPC=P
fi
echo "proc set hpc to $HPC"
NCPU=$(sed -n 's/^NCPU *= *\([0-9]*\).*/\1/p' $CONFIG)
if [[ -z $NCPU ]]; then
	NCPU=1
fi
MAXQS=$(sed -n 's/^MAXQS *= *\([0-9]*\).*/\1/p' $CONFIG)
if [[ -z $MAXQS ]]; then
	MAXQS=1
fi
date
time
# Setup BLAST database
echo "Setting up BLAST database ..."
if [[ ! -d $OID_USER_DIR/blastdb ]]; then
	echo "Sequence directory \"$OID_USER_DIR/blastdb\" does not exist ... exiting"
	exit 1
fi
$OID_BIN/setup_blastdb.sh "${INGROUP[@]}" "${OUTGROUP[@]}"
if [[ $? -ne 0 ]]; then
	echo "Unable to setup BLAST db ... exiting"
	exit 1
fi

date
time
# Create relevant directories
if [[ ! -d $OID_USER_DIR/blast ]]; then
	mkdir $OID_USER_DIR/blast
fi
if [[ ! -d $OID_USER_DIR/data ]]; then
	mkdir $OID_USER_DIR/data
fi
if [[ ! -d $OID_USER_DIR/log/job ]]; then
	mkdir -p $OID_USER_DIR/log/job
fi

# Generate job script
JOB_SCRIPT=$OID_USER_DIR/run_oid_job.sh
MKJOB=mkjob$HPC.sh
$OID_HOME/bin/$MKJOB $JOB_SCRIPT
# End job script
chmod a+x $JOB_SCRIPT

# All-all BLAST
echo hpc set $HPC
if [[ $SLURM_JOB_ID ]]; then
	echo proc job $SLURM_JOB_ID
	echo proc cpu $SLURM_CPUS_PER_TASK
	echo proc ntask $SLURM_NTASKS
	export SOFT=$SCRATCH/software/bigplant
fi
if [[ ! -s $OID_USER_DIR/blast/blastres.blst ]]; then
	cp $OID_USER_DIR/blastdb/combined.fa $OID_USER_DIR/blast
	#        source activate $SOFT
	#        module load perl
	$OID_HOME/bin/new_blast_parts.pl #make partn.faa (pgm estimates size
	#        NPART = $(/bina/ls OID_USER_DIR/blast/*.* | grep -c ".faa")
	$OID_HOME/bin/qsblast.pl -g 16 -n 12 -w 12 -q $MAXQS
	#fi
	#echo 'proc finished'
	if ! [[ -f $OID_USER_DIR/blast/blastres.blst && -f $OID_USER_DIR/blast/genelen.blst ]]; then
		echo "Error: failed to generate BLAST results database"
		exit 1
	else
		echo "Archiving part files..."
		tar -czf $OID_USER_DIR/blast/Parts.tar.gz $OID_USER_DIR/blast/Part[0-9]* --remove-files
		tar -czf $OID_USER_DIR/blast/speciesParts.tar.gz $OID_USER_DIR/blast/[A-z]*Part[0-9]* --remove-files
	fi
else
	echo "All-aginst-all BLAST results exist ... skipping"
fi
date
time
# Create gene families
#if ! /bin/ls -d $OID_USER_DIR/data/[1-9] >/dev/null 2>&1; then
if [[ ! -f $OID_USER_DIR/data/.family.done ]]; then
	echo "Creating families ..."
	#	if [[ -f $OID_USER_DIR/blast/clusters ]]; then
	# Clustering done, just create family directories
	#		$OID_HOME/bin/orthologid.pl -f
	$OID_HOME/bin/runmcl.pl $HIMEM 20 # $MY_MEM 20
	#	else
	#		JOBID=$(qsub -l nodes=1:ppn=$NCPU,walltime=12:00:00 $JOB_SCRIPT -v arg1="-f" | grep '^[0-9]')
	#		# Wait for clustering to finish
	#		while sleep 60; do
	#			if [[ $(qstat $JOBID 2>/dev/null | grep -c $JOBID) -eq 0 ]]; then
	#				break
	#			fi
	#		done
	#	fi
	if [[ $? -ne 0 ]]; then
		echo "Family clustering failed!"
		exit 1
	fi
else
	echo "Gene families already exist"
fi
echo "running new job select"
## new job to grep all data subdir for Family size
$OID_HOME/bin/rdfamdb.pl
#
rm log/job/schedone
while [[ ! -f log/job/schedone ]]; do
	$OID_HOME/bin/orthologid.pl -s 'hello'
	if [[ ! -f log/job/schedone ]]; then
		echo "orthologid.pl -s aborted before finish"
		date
	fi
done
date
time
# Extract orthologs
echo "Extracting orthologs ..."
$OID_BIN/orthologid.pl -O
date
time

#Fix tree run errors
echo "Checking for tree run errors..."
RERUN=0
source $OID_BIN/checktrees.sh
if [[ $RERUN -eq 1 ]]; then
	# Extract orthologs again
	echo "Extracting skipped orthologs ..."
	$OID_BIN/orthologid.pl -O
	date
	time
fi

# Generate big matrix
echo 'Generating matrix ...'
$OID_BIN/orth2matrix.pl

date
time

# Run tree searches
echo "Tree search ..."
if [[ -f jac.tre ]]; then
	echo "Tree file already exists"
else
	tnt bground p $OID_HOME/PostProcessing/mpt.proc
	tnt bground p $OID_HOME/PostProcessing/jac.proc
fi
while sleep 300; do
	#	if [[ ! -f Jacknife.tre ]]; then
	if [[ ! -f jac.tre ]]; then
		continue
	fi
	break
done
# Calculate PBS values
#echo "Calculating PBS values ..."
#$OID_BIN/pbs_split_slurm.pl

date
time

echo "Processing trees..."
python $OID_HOME/PostProcessing/fix_tree.py mpt.tre mpt_fixed.tre
python $OID_HOME/PostProcessing/fix_tree.py mpt.nel mpt_nel_fixed.tre
python $OID_HOME/PostProcessing/fix_tree.py jac.tre jac_fixed.tre

echo "Pipeline complete!"
date
time
