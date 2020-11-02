#!/usr/bin/perl -w
#
# translate2aa.pl
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
# Translate dna sequences in fasta format using EMBOSS getorf.
# The longest ORF is used.  Input sequences shorter than a
# threshold (minLen) is skipped.
#
# getorf must be located in PATH.
#
# Author: Ernest K. Lee <elee@amnh.org>
# $Id: translate2aa.pl 85 2008-05-27 05:01:16Z ernie $
#
use Bio::SeqIO;
use IPC::Open2;
use strict;

my $usage = "translate2aa.pl [ -m minlen ] dna_fasta_file";

my $minLen = 30;  # Minimum number of nucleotides to be reported

sub dna2aa($) {
        my $getorfProg = "getorf --filter -minsize $minLen";
        my $dna = shift @_;
        return if length($dna) < $minLen;

        # Run getorf
        my $pid = open2(*RFH, *WFH, $getorfProg);
        print WFH $dna;
        close WFH;
        my $maxLen = 0;
        my $aa;
        my $best = 0;
        while (<RFH>) {
                chomp;
                if (/^>.*\[(\d+)\D+(\d+)]/) {
                        my $len = abs($2-$1);
                        if ($len > $maxLen) {
                                $maxLen = $len;
                                $best = 1;
                                $aa = "";
                        }
                        else {
                                $best = 0;
                        }
                }
                elsif ($best) {
                        $aa .= $_;
                }
        }
        close RFH;
		waitpid($pid, 0);
        return $aa;
}

# Get options and arguments
if (@ARGV == 3 && $ARGV[0] eq "-m") {
	shift;
	$minLen = shift;
	die "Minimum length must be a positive integer!\n" if $minLen !~ /^\d+$/;
}
die "Usage: $usage\n" if @ARGV != 1;
my $fna = $ARGV[0];
die "$fna does not exist!\n" if ! -f $fna;

my $in = Bio::SeqIO->new(-file => $fna, -format => 'Fasta');
while (my $seq = $in->next_seq) {
        my $aa = dna2aa($seq->seq);
        if (!defined $aa) {
                warn "Could not translate ", $seq->id, "\n";
                next;
        }
        print ">", $seq->id, "\n";
        print "$aa\n";
}

