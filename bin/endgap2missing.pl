#!/usr/bin/perl -w
#
# endgap2missing.pl
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
# Changes leading and trailing gaps of sequences to '?'s.
#
# Author: Ernest Lee <elee@amnh.org>
#

use strict;

my @seqNames;
my %seqs;
while (<>) {
	chomp;
	if (/^>(\S+)/) {
		@seqNames = (@seqNames, $1);
		$seqs{$1} = "";
	}
	else {
		s/\s//g;
		$seqs{$seqNames[@seqNames-1]} .= $_;
	}
}

my $leading;
my $trailing;
foreach (@seqNames) {
	if ($seqs{$_} =~ /(^[\s-]*)([A-Z]+[A-Z-]*[A-Z])([\s-]*$)/i) {
		$leading = $1;
		$trailing = $3;
		$leading =~ tr/-/?/;
		$trailing =~ tr/-/?/;
		$seqs{$_} = $leading.$2.$trailing;
		print ">$_\n", $seqs{$_}, "\n";
	}
}

