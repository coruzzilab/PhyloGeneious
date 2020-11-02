#!/bin/bash
#
# Setup OrthologID BLAST db
#
# Copyright (C) 2006-2011 Ernest K. Lee
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
# $Id: setup_blastdb.sh 91 2009-03-02 06:24:06Z ernie $
#

echo "setup_blast+"
if [[ $OID_HOME == "" ]]; then
	echo -u2 "OID_HOME not defined ... exiting"
	exit 1
fi
if [[ $OID_USER_DIR == "" ]]; then
	echo -u2 "OID_USER_DIR not defined ... exiting"
	exit 1
fi
export PATH=~/blastplus/bin:$PATH

#FORMATDB=$(which formatdb 2>/dev/null)
#if [[ -z $FORMATDB ]]; then
#	BLASTDIR=$(grep BLASTHOME $OID_HOME/config | sed 's/.*=[[:space:]]*//')
#	if [[ -d $BLASTDIR ]]; then
#		echo -u2 "BLAST installation not found!  Make sure BLASTHOME is defined in config file or BLAST executables are in your PATH."
#		exit 1
#	else
#		FORMATDB=$BLASTDIR/bin/formatdb
#	fi
#fi
#echo "format $FORMATDB"
Combined="combined.fa"
cd $OID_USER_DIR/blastdb
if [[ -s $Combined ]]; then
    echo -u2 "combined fa already exists - just apending"
else
    >$Combined
fi

for i in "$@"; do
    echo $i
	
	# Check if BLAST DB exists
	if [[ -s $i.psq ]]; then
		echo "BLAST db file for $i already exists ... skipping"
		continue
	fi
	
	# Find fasta file
	fasta_file=""
	if [[ -f $i.seq ]]; then
		fasta_file=$i.seq
	elif [[ -f $i.faa ]]; then
		fasta_file=$i.faa
	elif [[ -f $i.fa ]]; then
		fasta_file=$i.fa
	elif [[ -f $i.fasta ]]; then
		fasta_file=$i.fasta
	elif [[ -f $i.pep ]]; then
		fasta_file=$i.pep
	else
		echo -u2 "Fasta file for $i does not exist"
		exit 1
	fi

	# Check deflines
	if grep -q '^>.*[^A-Za-z0-9#._]' $fasta_file; then
		echo -u2 "Error: Illegal characters in gene names in $fasta_file."
		echo -u2 "Only letter, digits, ., and _ are allowed, with # being the species delimiter."
		exit 1
	fi
	if grep -q '^>[^#]*$' $fasta_file; then
		echo -u2 "Warning: The following gene names do not contain a # character to delimit the species name and geneID:"
		grep '^>[^#]*$' $fasta_file
	fi
	
	# Create BLAST DB (assuming protein sequences)
	echo "Creating ${i%%.*} BLAST db"
    cat $fasta_file >>$Combined
# 	$FORMATDB -n $i -i $fasta_file -o T
    makeblastdb -dbtype 'prot' -in $fasta_file -title $i -out $i -parse_seqids
done

		
