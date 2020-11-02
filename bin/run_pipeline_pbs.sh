#!/bin/ksh
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
HIMEM="12GB"
##export OID_USER_DIR=`pwd`
##export OID_HOME=/home/cmz209/orthotnt/OID_nw3
if [[ -z $MY_SHELL ]]; then
	print -u2 "MY_SHELL not defined - won't restart properly"
    export MY_SHELL=$0
##	exit 1
fi
print "will restart $MY_SHELL"
$OID_HOME/bin/gettime.pl >starttime
#module load mercurial/intel/2.1.2
#module load mcl/intel/12-068
#module load mafft/intel/6.864

# Check environment variables and user directory
#if ! env | grep -q "^OID_HOME="; then
#	print -u2 "OID_HOME not defined ... exiting"
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
    print "copied procfiles.txt from $OID_HOME"
fi

OID_BIN=$OID_HOME/bin
print "oidbin $OID_BIN"

#if [[ $# -eq 0 ]] && ! env | grep -q "^OID_USER_DIR="; then
#	print -u2 "Must specify run directory as argument or define OID_USER_DIR"
#	exit 1
#elif [[ $# -eq 1 ]]; then
#	export OID_USER_DIR="$1"
#elif [[ $# -gt 1 ]]; then
#	print -u2 "Usage: $0 [ OID_USER_DIR ]"
#	exit 2
#fi
#if [[ ! -d $OID_USER_DIR ]]; then
#	print -u2 "OID_USER_DIR \"$OID_USER_DIR\" is not a directory"
#	exit 1
#fi

# Use dir name as run name
OID_RUN=$(basename $OID_USER_DIR)

print "== Starting OrthologID version $OID_VERSION pipeline =="
print "Run directory is $OID_USER_DIR"

# Check/parse config file
CONFIG=$OID_USER_DIR/config
if [[ ! -f $CONFIG ]]; then
    print -u2 "Config file required -- a sample is $OID_HOME/config/config "
	print -u2 "... exiting"
	exit 1
fi
set -A INGROUP $(grep INGROUP $CONFIG | cut -f2 -d=)
set -A OUTGROUP $(grep OUTGROUP $CONFIG | cut -f2 -d=)
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
print "Setting up BLAST database ..."
if [[ ! -d $OID_USER_DIR/blastdb ]]; then
	print -u2 "Sequence directory \"$OID_USER_DIR/blastdb\" does not exist ... exiting"
	exit 1
fi
$OID_BIN/setup_blastdb.sh "${INGROUP[@]}" "${OUTGROUP[@]}"
if [[ $? -ne 0 ]]; then
	print -u2 "Unable to setup BLAST db ... exiting"
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
cat <<EOF >$JOB_SCRIPT
#!/bin/sh
#
# PBS job script for orthologid
#

#PBS -S /bin/bash
#PBS -j oe
#PBS -l mem=$HIMEM
##PBS -o /scratch/cmz209/orthotnt/oidTest9/log/job/
#PBS -o log/job/
#PBS -q $PBSQ
#PBS -N $OID_RUN
#PBS -V 

PATH=$OID_BIN:$PATH

cd \$OID_USER_DIR
date
time
orthologid.pl "\$arg1" "\$arg2"

date
time
if [[ "\$arg1" == "-b" ]]; then
	touch blast/.\$arg2.done
fi
EOF
# End job script
chmod a+x $JOB_SCRIPT

# All-all BLAST
if [[ ! -s $OID_USER_DIR/blast/blastres.blst ]]; then
        cp $OID_USER_DIR/blastdb/combined.fa $OID_USER_DIR/blast
        $OID_HOME/bin/new_blast_parts.pl    #make partn.faa (pgm estimates size
#        NPART = $(/bina/ls OID_USER_DIR/blast/*.* | grep -c ".faa")
        $OID_HOME/bin/qsblast.pl  -g 16 -n 12 -w 12 -q $MAXQS
#	for i in "${INGROUP[@]}" "${OUTGROUP[@]}" ; do
##		print "Submitting job for $i-against-all BLAST..."
#		qsub -l nodes=1:ppn=$NCPU,walltime=8:00:00 $JOB_SCRIPT -v arg1="-b",arg2="$i"
#	done
#	print "Waiting for BLAST jobs to finish ..."
#date
#time
#	while sleep 300; do
#		for i in "${INGROUP[@]}" "${OUTGROUP[@]}" ; do
#			if [[ ! -f $OID_USER_DIR/blast/.$i.done ]]; then
#				continue 2
#			fi
#		done
#		break
#	done
#date
#time

#	print "Merging BLAST data ..."
#	$OID_BIN/orthologid.pl -B
	if ! [[ -f $OID_USER_DIR/blast/blastres.blst && -f $OID_USER_DIR/blast/genelen.blst ]]; then
		print -u2 "Error: failed to generate BLAST results database"
		exit 1
	fi
else
		print "All-aginst-all BLAST results exist ... skipping"
fi
date
time
# Create gene families
if ! /bin/ls -d $OID_USER_DIR/data/[1-9] >/dev/null 2>&1; then
	print "Creating families ..."
#	if [[ -f $OID_USER_DIR/blast/clusters ]]; then
		# Clustering done, just create family directories
		$OID_HOME/bin/orthologid.pl -f
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
		print -u2 "Family clustering failed!"
		exit 1
	fi
else
	print "Gene families already exist"
fi
print "running new job select"
## new job to grep all data subdir for Family size
$OID_HOME/bin/rdfamdb.pl
#
rm log/job/schedone
while [[ ! -f log/job/schedone ]]; do
$OID_HOME/bin/orthologid.pl -s 'hello'
if [[ ! -f log/job/schedone ]]; then
   print "orthologid.pl -s aborted before finish"
   date
fi
done;
date
time
# Extract orthologs
print "Extracting orthologs ..."
$OID_BIN/orthologid.pl -O
date
time

# Generate big matrix
print 'Generating matrix ...'
$OID_BIN/orth2matrix.pl

date
time
