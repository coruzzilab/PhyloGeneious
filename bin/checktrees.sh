#!/bin/bash
#
#
# Script for running OrthologID pipeline
# This script makes use of PBS for job submission.
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
# Author: Veronica M. Sondervan

runtreejob () {
if [[ $HPC -eq "S" ]]; then
	sbatch --mem=8GB --time=12:00:00 --wrap="tnt 'p oid.proc;zzz'"
else
	echo "tnt 'p oid.proc;zzz'" | qsub -l mem=8GB,walltime=12:00:00 -N wrap
fi
}

checkqueue () {
if [[ $HPC -eq "S" ]]; then
	squeue -u $USER | grep -c wrap
else
	qstat -u $MYUSER | grep -c wrap
fi
}

JOBFAMS=()

TREEPROG=$(sed -n 's/^TREEPROGRAM *= *\([TO].*\)/\1/p' $OID_USER_DIR/config)
if [[ -z $TREEPROG ]]; then
	TREEPROG="TNT"
fi

if [[ $TREEPROG -eq "TNT" ]]; then
	for FAMILYDIR in `comm -23 <(ls $OID_USER_DIR/data/*/FAMILY | grep -v "S" | sed "s/FAMILY//g"| sort) <(ls $OID_USER_DIR/data/*/oid.tre | grep -v "S" | sed "s/oid.tre//g" | sort)`; do #get fams with missing trees
		if [[ $(grep "Increase slack" ${FAMILYDIR}*log | tail -n1|wc -c) -gt 0 ]]; then
			echo "re-attempting tree for ${FAMILYDIR}"
			RERUN=1
			echo "re-attempting tree for ${FAMILYDIR}" > $OID_USER_DIR/log/treefix.log
			cd $FAMILYDIR
			while [ ! -f ${FAMILYDIR}oid.tre ]; do
				SLACKVAL=$(($(grep "Increase slack" ${FAMILYDIR}oid.log|tail -n1|cut -d" " -f6)+10))
				sed -i "s/slack [0-9]*/slack ${SLACKVAL}/g" ${FAMILYDIR}oid.proc
				tnt 'silent=console;silent=buffer;p oid.proc;zzz' &> $OID_USER_DIR/log/treefix.log
				echo ""
			done
			echo "built tree for ${FAMILYDIR}"
		elif [[ $( grep "out of ram" ${FAMILYDIR}*log | tail -n1|wc -c) -gt 0 ]]; then
			echo "Family in ${FAMILYDIR} is too large: skipping for orthology"
		#	cd $FAMILYDIR
		#	JOBFAMS+=(${FAMILYDIR})
		#	while [ ! -f ${FAMILYDIR}oid.tre ]; do
				#edit when figure out largest family run fix
				#LOG=`runbigtreejob`
				#echo $FAMILYDIR $LOG
		#	done
		else
			echo "Unknown tree building error for ${FAMILYDIR}: will attempt to rebuild"
			cd $FAMILYDIR
			if [[ $(grep -c "drift" ${FAMILYDIR}oid.proc) -gt 0 ]]; then
				JOBFAMS+=(${FAMILYDIR})
				LOG=`runtreejob`
				echo $FAMILYDIR $LOG
			else
				echo "re-attempting tree for ${FAMILYDIR}" > $OID_USER_DIR/log/treefix.log
				tnt 'silent=console;silent=buffer;p oid.proc;zzz' &> $OID_USER_DIR/log/treefix.log
				echo ""
				if [ ! -f ${FAMILYDIR}oid.tre ]; then
					echo "Rebuild for ${FAMILYDIR} failed: skipping for orthology"
				else
					RERUN=1
					echo "built tree for ${FAMILYDIR}"
				fi
			fi
		fi
	done
fi

cd $OID_USER_DIR
if [[ ${#JOBFAMS[@]} -gt 0 ]]; then
	NJOBS=`checkqueue`
	while [[ $NJOBS -gt 0 ]]; do
		sleep 300
		NJOBS=`checkqueue`
	done
	for FAMILYDIR in ${JOBFAMS[@]}; do
		if [ ! -f ${FAMILYDIR}oid.tre ]; then
			echo "Rebuild for ${FAMILYDIR} failed: skipping for orthology"
		else
			RERUN=1
			echo "built tree for ${FAMILYDIR}"
		fi
	done
fi
