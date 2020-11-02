#!/usr/bin/perl -w
#
# bigMatrix.pl
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
# Generate supermatrix, with partitions
# by number of species present, GO terms and EC numbers
# Output is bigmatrix.nex.
#
# Author: Ernest K Lee <elee@amnh.org>
#

my $OID_HOME;

BEGIN {
	$OID_HOME = $ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n"
	  if !defined($OID_HOME);
}

use lib "$OID_HOME/lib";
use OrthologID;
use OrthologID::Utils;
use Cwd;
use strict;

my $OID_DATADIR     = "$OID_HOME/data";
my $matrixFile      = "bigmatrix.nex";    # output matrix file
my $ecFile          = "ec.list";          # gene-EC# table
my $goFile          = "go.list";          # gene-GO terms table
my $ecMapFile       = "ec_map.tab";       # EC#-EC map# table (KEGG)
my $pathwayNameFile = "map_title.tab";    # EC map#-pathway table (KEGG)

my $savDir = getcwd;

my @INGROUP  = getIngroup();
my @OUTGROUP = getOutgroup();
my %taxa;    # taxa to include in matrix (1=include; 0=exclude)
foreach ( @INGROUP, @OUTGROUP ) {
	$taxa{$_} = 1;
}
my @ogroup  = ();    # list of hashes of alignments
my @spCount = ();    # species counts for each ortholog group (partition)

# Read in EC and GO tables
my %ec;              # keys are EC#s; vals are pointers to list of OG#s
my @unassignedEC;    # list of OG#s without EC assignments
open EFH, $ecFile or die "Failed to open EC table: $!";
my @ecList = <EFH>;
chomp(@ecList);
close EFH;

my %go;              # keys are GO terms; vals are pointers to list of OG#s
my @unassignedGO;    # list of OG#s without GO assignments
open GFH, $goFile or die "Failed to open GO table: $!";
my @goList = <GFH>;
chomp(@goList);
close GFH;

my @noAnno;          # partitions with no annotation (no GO nor EC)

# Read ec-pathway map
my %path2ec;         # pathway num - EC num list (ref) map
open MFH, $ecMapFile or die "Failed to open EC map: $!";
while (<MFH>) {
	chomp;
	/\t/;
	my $ec = $`;
	my @paths = split( ' ', $' );
	foreach my $path (@paths) {
		my $ecl = $path2ec{$path};    # list of ec numbers for "this" pathway
		if ($ecl) {
			push( @$ecl, $ec );
		}
		else {
			$path2ec{$path} = [ ($ec) ];
		}
	}
}
close MFH;

# Read pathway-name map
my %pathwayNum;                       # pathway number - pathway name map
open TFH, $pathwayNameFile or die "Failed to open map title file: $!";
while (<TFH>) {
	chomp;
	/\t/;
	$pathwayNum{$'} = $`;
}
close TFH;

chdir $OID_DATADIR;
foreach my $dir (<[1-9]*>) {
	open OFH, "$dir/orthologs" or next;
	my $seqH = readFasta("$dir/FAMILY.aligned");
	print "Processing family $dir ...\n";

	# if only one species has more than one gene present, then include 'em all
	# this is to avoid breaking up of probable ortholog sets
	# should be fixed to make it more robust
	my %t;
	my $multiGeneSp = 0;
	foreach my $gene (keys %$seqH) {
		$gene =~ /#/;
		$t{$`}++;
		$multiGeneSp++ if $t{$`} > 1;
	}
	if ( $multiGeneSp == 1) {
		print " ** include all **\n";
		my $tmpListFile = "/tmp/genmat.tmplist$$";
		close OFH;
		open OFH, ">$tmpListFile" or die "Cannot create tmp file: $!\n";
		my @tlist = keys %$seqH;
		print OFH "@tlist\n";
		close OFH;

		# reopen for reading
		open OFH, $tmpListFile or die "Cannot reopen tmp file: $!\n";
	}
	while (<OFH>) {
		chomp;
		my @olist = split;
		my @taxlist = map { /^[^#]*/; $& } @olist;
		print " set = { @taxlist }\n";
		my %numSp = ();
		my %seq;
		foreach my $gene (@olist) {

			# Only care about species, not strains
			#$gene =~ /^(.{4}).*#(.*)/;
			$gene =~ /#/;
			my $sp     = $`;
			my $geneID = $';
			$seq{$gene} = $$seqH{$gene};

			# $seq{$sp."#".$geneID} = $$seqH{$gene};
			$numSp{$sp}++ if $taxa{$sp} == 1;
		}
		if ( keys(%numSp) > 1 ) {

			# store partition
			push( @ogroup, \%seq );

			# record species count
			push( @spCount, scalar( keys(%numSp) ) );

			my $partID = "OG" . scalar(@ogroup);
			my $noEC = 1;    # keep track of ortholog groups without EC/GO
			my $noGO = 1;
			my %ECcnt;       # weight of genes annotated with specific EC's
			my %GOcnt;       # weight of genes annotated with specific GO's
			my $totAnnoWgt = 0;   # total possible weight in this ortholog group
			foreach my $gene (@olist) {
				my $line;

				# Get EC numbers
				($line) = grep( /^$gene/, @ecList );
				$line = "" if !defined($line);
				my @ecl = split( ' ', $line );
				shift @ecl;
				print "  EC: @ecl <- $gene\n" if @ecl > 0;
			  EC: foreach (@ecl) {					
					my @newList = ();
					if ( defined( $ec{$_} ) ) {
						@newList = @{ $ec{$_} };
					}
					foreach my $n (@newList) {

						# ortholog group already assigned this EC
						next EC if $n eq $partID;
					}
					push( @newList, "$partID" );
					$ec{$_} = \@newList;
					$noEC = undef;
				}
			}
			push( @unassignedEC, $partID ) if $noEC;

			foreach my $gene (@olist) {
				my $line;

				# Get GO terms
				($line) = grep( /^$gene/, @goList );
				$line = "" if !defined($line);
				my @gol = split( ' ', $line );
				shift @gol;
				print "  GO: @gol <- $gene\n" if @gol > 0;
			  GO: foreach (@gol) {
					my @newList = ();
					if ( defined( $go{$_} ) ) {
						@newList = @{ $go{$_} };
					}
					foreach my $n (@newList) {

						# ortholog group already assigned this GO
						next GO if $n eq $partID;
					}
					push( @newList, "$partID" );
					$go{$_} = \@newList;
					$noGO = undef;
				}
			}
			push( @unassignedGO, $partID ) if $noGO;

			# Keep track of no EC *and* no GO
			push( @noAnno, $partID ) if $noEC && $noGO;

		}
	}

}

chdir $savDir;
print "Generating matrix from " . scalar(@ogroup) . " ortholog groups ...\n";
my @taxonList;
while ( my ( $taxon, $include ) = each %taxa ) {
	push( @taxonList, $taxon ) if $include == 1;
}
my $matrix = genSuperMatrix( [ sort @taxonList ], @ogroup );
open FH, ">$matrixFile" or die "Cannot write matrix file: $!";
print FH $matrix;

# Write charsets by species present
print "Writing charsets by number of species ...\n";
my %partBySp;
foreach ( 0 .. $#ogroup ) {
	$partBySp{ $spCount[$_] } .= "OG" . ( $_ + 1 ) . " ";    # OG# starts with 1
}
print FH "[ Partitions by number of species in ortholog group ]\n";
print FH "BEGIN SETS;\n\n";
foreach ( sort { $a <=> $b } keys %partBySp ) {
	print FH "CHARSET SP$_ = $partBySp{$_};\n";
}
print FH "END;\n\n";

print FH "BEGIN SETS;\n\n";

# Write charsets by EC numbers
print "Writing charsets by EC numbers ...\n";
print FH "[ Partitions by EC number ]\n";
foreach ( sort keys %ec ) {
	my $ecMod = $_;
	$ecMod =~ s/-/_/g;
	print FH "CHARSET E$ecMod = @{$ec{$_}};\n";
}
print FH "CHARSET Enone = @unassignedEC;\n" if @unassignedEC > 0;
print FH "\n\n";

# Write charsets by pathways
print "Writing charsets by pathways ...\n";
print FH "[ Partitions by pathways ]\n";
foreach ( sort keys %pathwayNum ) {
	my $pNum  = $pathwayNum{$_};    # KEGG pathway#
	my $pName = $_;
	$pName =~ tr$ ,'()/-$_$s;       # replace illegal characters (for PAUP)'
	my %ogl;
	foreach my $e ( @{ $path2ec{$pNum} } ) {    # list of EC#s for this pathway
		if ( $ec{$e} ) {
			foreach my $og ( @{ $ec{$e} } ) {    # list of OG#s for this EC#
				$ogl{$og} = 1;
			}
		}
	}
	if ( keys %ogl > 0 ) {
		my @ogl = keys %ogl;
		print FH "CHARSET $pName = @ogl;\n";
	}
}
print FH "\n\n";

# Write charsets by GO terms
print "Writing charsets by GO terms ...\n";
print FH "[ Partitions by GO term ]\n";
foreach ( sort { $a <=> $b } keys %go ) {
	print FH "CHARSET G$_ = @{$go{$_}};\n";
}
print FH "CHARSET Gnone = @unassignedGO;\n" if @unassignedGO > 0;
print FH "\n\n";

print FH "CHARSET EGnone = @noAnno;\n" if @noAnno > 0;
print FH "\n\n";

print FH "END;\n\n";

close FH;

