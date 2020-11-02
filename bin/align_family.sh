#!/bin/bash
#
# align_family.sh
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
# Align using MAFFT with default options.
# Must be run in top level of a gene family directory.
#
# Author: Ernest Lee <elee@amnh.org>
#

if [[ $# -ne 2 ]] ; then
	echo "Usage: align_family fasta-file output-file"
	exit 2
fi

INPUT=$1
OUTPUT=$2

# Check I/O files and directory
if [[ ! -r $INPUT ]] ; then
	echo -u2 "$INPUT does not exist or is not readable"
	exit 1
fi
if [[ -s $OUTPUT ]] ; then
	echo -u2 "$OUTPUT exists"
	exit 1
fi

# Check MAFFT
if ! which mafft >/dev/null 2>/dev/null; then
	echo -u2 "MAFFT is not installed or not found."
	exit 1
fi

# MAFFT alignment
numSeq=$(grep -c '^>' $INPUT)
if ((numSeq<500)) ; then
	MAFFT_EXE='mafft-linsi'
else
	# Let MAFFT decides for large family
	MAFFT_EXE='mafft --auto'
fi

echo -n "Using $MAFFT_EXE ... "
$MAFFT_EXE --anysymbol --quiet $INPUT > $OUTPUT
status=$?

# Re-run mafft if output is empty...
if [[ ! -s $OUTPUT ]]; then
	echo -n "failed\nRe-running mafft ... "
	MAFFT_EXE='mafft --auto'
	$MAFFT_EXE --anysymbol --quiet $INPUT > $OUTPUT
	status=$?
fi

echo "done"

exit $status

