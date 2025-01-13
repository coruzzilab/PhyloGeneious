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

if [[ $SEARCHTYPE -eq 'B' ]]; then
    echo "setup_blast+"
elif [[ $SEARCHTYPE -eq 'D' ]]; then
    echo "setup_diamond"
elif [[ $SEARCHTYPE -eq 'M' ]]; then
    echo "setup mmseqs2"
# elif [[ $SEARCHTYPE -eq 'K' ]]; then
#     echo "setup kaamer"
else
    echo "search method not recognized. please use blastp, diamond, or mmseqs2."
    exit 1
fi
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
	if [[ $SEARCHTYPE -eq 'B' ]]; then
            if [[ -s $i.psq ]]; then
		echo "BLAST db file for $i already exists ... skipping"
		continue
	    fi
	elif [[ $SEARCHTYPE -eq 'D' ]]; then
    	    if [[ -s $i.dmnd ]]; then
	        echo "DIAMOND db file for $i already exists ... skipping"
		continue
	    fi
	elif [[ $SEARCHTYPE -eq 'M' ]]; then
            if [[ -s $i.DB ]]; then
	        echo "MMseqs2 db file for $i already exists ... skipping"
		continue
	    fi
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
	if grep -v '^>' $fasta_file | grep -q '[^AaBbCcDdEeFfGgHhIiKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*]'; then
 		echo -u2 "Error: Illegal characters detected in protein sequences in $fasta_file."
   		grep -v '>' $fasta_file | grep '[^AaBbCcDdEeFfGgHhIiKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz*]'
     		exit 1
       fi

	# Create BLAST DB (assuming protein sequences)
	if [[ $SEARCHTYPE -eq 'B' ]]; then
    	    echo "Creating ${i%%.*} BLAST db"
	    cat $fasta_file >>$Combined
	#     $FORMATDB -n $i -i $fasta_file -o T
	    makeblastdb -dbtype 'prot' -in $fasta_file -title $i -out $i -parse_seqids
	elif [[ $SEARCHTYPE -eq 'D' ]]; then
            echo "Creating ${i%%.*} DIAMOND db"
	    cat $fasta_file >>$Combined
	    diamond makedb --in $fasta_file --db $i
	elif [[ $SEARCHTYPE -eq 'M' ]]; then
    	    echo "Creating ${i%%.*} MMseqs2 db"
	    cat $fasta_file >>$Combined
	    mmseqs createdb $fasta_file $i.DB --dbtype 1
	    mmseqs/bin/mmseqs createindex $i.DB tmp
	# elif [[ $SEARCHTYPE -eq 'K' ]]; then
    	#     echo "Creating ${i%%.*} kaamer db"
	#     cat $fasta_file >>$Combined
	#     go/bin/kaamer-db -make -f faa -i $fasta_file -d $i
	fi
done
