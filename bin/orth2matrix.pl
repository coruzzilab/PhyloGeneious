#!/usr/bin/env perl
#
# orth2matrix.pl
#
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
# Generate matrix from OrthologID orthologs files.
#
# Authors: Ernest K Lee <elee@amnh.org>
# Authors: Kranthi K Varala	<kvarala@purdue.edu>
#

use warnings;

my ( $OID_HOME, $OID_USER_DIR );

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);

}

use Cwd;
use Getopt::Std;
use lib "$OID_HOME/lib";
use OrthologID;
use OrthologID::Utils;
use strict;

# defaults
my $minNonMissingTaxa = 2;
my %excludeTaxon;    # Exclude taxon (key) if defined
my $needOutgroup
  ;    # At least one outgroup must be present in each partition if defined
my $MatrixFile             = "$OID_USER_DIR/Matrix.nex";
my $TNTMatrixFile          = "$OID_USER_DIR/Matrix.tnt";
my $LikelihoodMatrixFile   = "$OID_USER_DIR/Matrix.phy";
my $PartitionFile          = "$OID_USER_DIR/partitions.txt";
my $PartitionMembers       = "$OID_USER_DIR/partitionMembers.txt";
my $PartitionMembersUnfilt = "$OID_USER_DIR/partitionMembers_paralogs.txt";
my $alignedFile            = "FAMILY.aligned";
my %PhylipMat;    # Holds the phylip matrix for RAxML
my $pnames;       #List of partition names
$Getopt::Std::STANDARD_HELP_VERSION = 1;

sub HELP_MESSAGE() {
    print
"Usage: orth2matrix.v2.pl [ -O ] [ -m min_non_missing_taxa ] [ -x comma_delimited_taxon_list_to_exclude ]\n";
    print
"Setting the -O flag will remove all partitions that do not contain a sequence from the outgroup. Default: off\n";
    print
"-m sets the minimum number of taxa that need to be represented in a partition to include it. Default: 4\n";
    print
      "-x may be used to exclude taxa from the final matrix. Default: None\n";
}

# Get options
our ( $opt_m, $opt_O, $opt_x );
getopts('Om:x:');
$minNonMissingTaxa = $opt_m if $opt_m;
$needOutgroup      = 1      if $opt_O;
if ($opt_x) {
    foreach ( split( ',', $opt_x ) ) {
        $excludeTaxon{$_} = 1;
        print "Excluding $_ from matrix\n";
    }
}

# Check if matrix.nex already exists
die "\"$MatrixFile\" already exists ... will not overwrite\n" if -f $MatrixFile;
die "\"$TNTMatrixFile\" already exists ... will not overwrite\n"
  if -f $TNTMatrixFile;
die "\"$LikelihoodMatrixFile\" already exists ... will not overwrite\n"
  if -f $LikelihoodMatrixFile;

# Get taxa
my @INGROUP  = getIngroup();
my @OUTGROUP = getOutgroup();
my %isOutgroup;
my @incINGROUP;     # ingroup not excluded
my @incOUTGROUP;    # outgroup not excluded
my %taxa;           # taxa to include in matrix if defined
foreach (@INGROUP) {

    if ( !$excludeTaxon{$_} ) {
        $taxa{$_} = 1;
        push @incINGROUP, $_;
    }
}
foreach (@OUTGROUP) {
    $isOutgroup{$_} = 1;
    if ( !$excludeTaxon{$_} ) {
        $taxa{$_} = 1;
        push @incOUTGROUP, $_;
    }
}
my @ogroup = ();    # list of hashes of alignments
my @oanno  = ();    # list of annotations

my $savDir      = getcwd;
my $OID_DATADIR = "$OID_USER_DIR/data";
chdir $OID_DATADIR;
my $ogNum = 1;
my @fams = sort { $a <=> $b } <[0-9]*>;
foreach my $dir (@fams) {    #(<[1-9]*>) {

    my $hasOutgroup = 0;
    next if !-r "$dir/orthologs";
    print "Processing family $dir\n";
    my $seqH = readFasta("$dir/$alignedFile");
    die "Failed to read/parse alignment.\n" if !defined $seqH;

    # Determine if there is a global set of outgroup orthologs
    open OFH, "$dir/orthologs"
      or die "Cannot open orthologs file in \"$dir\": $!\n";
    my @globalOutgroup    = ();
    my %globalOutgroupSeq = ();
    my $totGene           = 0;
    while (<OFH>) {
        chomp;
        my @olist = split;
        my @taxlist = map { /^[^#]*/; $& } @olist;
        $totGene += @olist;
        my $allOutgroup = 1;    # only outgroup in this set
        my $allIngroup  = 1;    # only ingroup in this set
        my %seq;
        foreach my $gene (@olist) {

            $gene =~ /#/;
            my $sp = $`;
            $seq{$gene} = $$seqH{$gene};

            if ( $isOutgroup{$sp} ) {
                $allIngroup = 0;
            }
            else {
                $allOutgroup = 0;
            }
        }
        if ( !$allIngroup && !$allOutgroup ) {

            # do not use global outgroup
            @globalOutgroup = ();
            last;
        }
        if ( $allOutgroup && @olist > @globalOutgroup ) {  # pick the bigger set
            @globalOutgroup    = @olist;
            %globalOutgroupSeq = %seq;
        }
    }
    close OFH;
    if ( @globalOutgroup && $totGene == keys(%$seqH) ) {

        #print " global outgroup = { @globalOutgroup }\n";
    }
    else {
        @globalOutgroup = ();
    }

    open OFH, "$dir/orthologs"
      or die "Cannot open orthologs file in \"$dir\": $!\n";
    while (<OFH>) {
        chomp;
        my @olist = split;
        my @taxlist = map { /^[^#]*/; $& } @olist;
        print " set = { @taxlist }\n";
        my %numSp = ();
        my %seq;
        my $numIngroup  = 0;    # number of ingroup taxa
        my $allOutgroup = 1;    # only outgroup in this set

        foreach my $gene (@olist) {

            $gene =~ /#/;
            my $sp     = $`;
            my $geneID = $';
            $seq{$gene} = $$seqH{$gene};

            # $seq{$sp."#".$geneID} = $$seqH{$gene};
            $numSp{$sp}++ if $taxa{$sp};
            if ( $numSp{$sp} == 1 && !$isOutgroup{$sp} ) {

                # first encounter of this ingroup taxon - increment
                $numIngroup++;
            }

            if ( $isOutgroup{$sp} ) {
                $hasOutgroup = 1;
            }
            else {
                $allOutgroup = 0;
            }
        }

        # count no. species in global outgroup
        my $numSpGlobalOutgroup = 0;
        my %numSpGOG;
        foreach my $gene (@globalOutgroup) {
            $gene =~ /#/;
            $numSpGOG{$`}++;
        }

        #		if ($allOutgroup && @globalOutgroup > 0 ||
        #			keys(%numSp) + keys(%numSpGOG) < $minNonMissingTaxa ||
        #			$needOutgroup && ! $hasOutgroup) {
        if (   $numIngroup < $minNonMissingTaxa
            || $needOutgroup && !$hasOutgroup )
        {
            print " ...skipped\n";
        }
        else {
            # store partition and annotation
            if (@globalOutgroup) {
                print " +global outgroup\n";
                while ( my ( $g, $s ) = each %globalOutgroupSeq ) {
                    $seq{$g} = $s;
                }
                push( @olist, @globalOutgroup );
            }
            push( @ogroup, \%seq );
            push( @oanno,  "OG$ogNum|Family$dir|@olist" );
            print "OG$ogNum: @olist\n";
            $ogNum++;
        }
    }

    close OFH;
}

chdir $savDir;

print "Generating matrix from " . scalar(@ogroup) . " ortholog groups ...\n";

# Sort taxa and generate matrix
my @taxaSorted = sort keys %taxa;
my $matrix = genSuperMatrix( \@taxaSorted, \@incOUTGROUP, \@oanno, @ogroup );

open MFH, ">$MatrixFile"
  or die "Cannot open file for writing PAUP matrix: $!\n";
print MFH "$matrix\n";

# Add partition size charsets

# Define outgroup
if ( @incOUTGROUP > 0 ) {
    print MFH "BEGIN PAUP;\n";
    print MFH "outgroup @incOUTGROUP;\n";
    print MFH "END;\n";
}

close MFH;

open TNT, ">$TNTMatrixFile"
  or die "Cannot open file for writing TNT matrix: $!\n";
open PHY, ">$LikelihoodMatrixFile"
  or die "Cannot open file for writing Likelihood matrix: $!\n";
open PART, ">$PartitionFile"
  or die "Cannot open file for writing partition length information: $!\n";
open PSET, ">$PartitionMembers"
  or die "Cannot open file for writing partition membership information: $!\n";
open PFSET, ">$PartitionMembersUnfilt"
  or die
  "Cannot open file for writing full partition membership information: $!\n";

print TNT "nstates 32;\ntaxonomy=;\nsilent = buffer ;\nxread\n";
my @mat = split /\n/, $matrix;
my ( $NTAX, $NCHAR );

# Parse the nexus header to retrieve number of taxa and number of characters
for ( my $i = 0 ; $i <= 5 ; $i++ ) {
    if ( $mat[$i] =~ /^DIMENSIONS/ ) {
        $mat[$i] =~ /NTAX=(\d+)\sNCHAR=(\d+)/;
        $NTAX  = $1;
        $NCHAR = $2;
    }
}

#print TNT "$NCHAR $NTAX\n";
print PHY " $NTAX $NCHAR\n";
open LOG, ">MatrixRecoding.log";

my $characters = 0;
my @curOG;
my @recodedOG;
my $part = 1;
my $pMembers;
my $pfMembers;
for ( my $i = 6 ; $i <= $#mat ; $i++ ) {
    if ( $mat[$i] =~ /^;/ ) { last; }
    if ( $mat[$i] =~ /^$/ ) { next; }
    if ( $mat[$i] =~ /^\[/ ) {
        if ( $#curOG > 1 ) {

            #print TNT "&[prot nogaps]\n";
            #print TNT @curOG;
            my $collapsedOG = join "", @curOG;
            push @recodedOG, code( $collapsedOG, $pMembers );

            #print TNT @recodedOG;
            $pMembers  = "";
            $pfMembers = "";
            undef @curOG;
        }
        my @OG = split /\s/, $mat[$i];
        my @GeneSet = ( split /\s/, ( split /\|/, $mat[$i] )[2] );
        my %oSet;
        foreach my $gene (@GeneSet) {
            my ( $taxon, $geneID ) = split /#/, $gene;
            $pfMembers .= "$taxon#$geneID ";
            unless ( defined $oSet{$taxon} ) {
                $pMembers .= "$taxon#$geneID ";
                $oSet{$taxon} = $geneID;
            }
        }
        my $partNumber = ( split /\|/, $OG[1] )[0];
        my $FamNumber  = ( split /\|/, $OG[1] )[1];
        $partNumber =~ s/OG/p/;
        my $partLength = $OG[-2];
        print PART "WAG, $partNumber = $partLength\n";
        print PSET $partNumber,  " = $pMembers\t$FamNumber\n";
        print PFSET $partNumber, " = $pfMembers\t$FamNumber\n";
        next;
    }
    my @curLine = split /\s+/, $mat[$i];
    $PhylipMat{ $curLine[0] } .= "$curLine[1] ";

    #if($isOutgroup{$curLine[0]}){$curLine[0]="\@OUTGROUP_".$curLine[0];}
    #else{$curLine[0]="\@INGROUP_".$curLine[0];}
    push @curOG, join "      ", @curLine, "\n";
}

#print TNT "&[prot nogaps]\n";
#print TNT @curOG;
my $collapsedOG = join "", @curOG;
push @recodedOG, code( $collapsedOG, $pMembers );
print TNT "$characters $NTAX\n";
print TNT @recodedOG;
print TNT ";\ncc - .;\n";
print TNT "outgroup[OUTGROUP;\n\n";
print TNT "proc/;";
close TNT;
close LOG;

foreach my $sps ( sort keys %PhylipMat ) {
    print PHY "$sps\t$PhylipMat{$sps}\n";
}
close PHY;
close PSET;
close PFSET;
close PART;

sub code {
    my @lines     = split( /\n/, $_[0] );
    my $partition = $_[1];
    my $r         = 0;                      #
    my $maxStates = 32;
    ############################### MAKE MATRIX
#if($#lines != ($NTAX-1)){
#	print(LOG "The number of taxa $#lines differs from expected partitions ($NTAX)! Partition '$partition' has been omitted.\n");
#	return();
#}
    my $matrix = ();
    $matrix->[0][0] = ();
    my $chars = ();
    for ( my $k = $#lines ; $k >= 0 ; $k-- ) {
        $lines[$k] =~ tr/[A-Z][a-z][0-9]\_?\- \*//cd
          ; #Remove any character that is not a letter, number, space or one of: _,?,-,*
        $lines[$k] =~ tr/ / /s;    # Replace multiple spaces with one space.
        my @line = split( / /, $lines[$k] );
        if ( length( $line[0] ) && length( $line[1] ) ) {
            if ( keys(%isOutgroup) ) {
                if ( exists( $isOutgroup{ $line[0] } ) ) {
                    $matrix->[$k][0] = '@OUTGROUP_' . $line[0];
                }
                else {
                    $matrix->[$k][0] = '@INGROUP_' . $line[0];
                }
            }
            else {
                $matrix->[$k][0] = $line[0];
            }

            #$matrix->[$k][0] = $line[0];
            my $seq = uc( join( '', @line[ 1 .. $#line ] ) );
            $seq =~ tr/\*X/\?\?/;    # Replace * and X (stop codons) with ?
	    $seq =~ tr/B/D/;    # Replace B (D or N) with D (most common)
	    $seq =~ tr/Z/E/;    # Replace Z (E or Q) with E (most common)
            $seq =~ tr/ACDEFGHIKLMNOPQRSTUVWY\?\-//cd
              ;    #Allow only AA letters or ? or - in alignment
            my @aa = split( //, $seq );
            if ( $chars && ( $chars != $#aa ) ) {
                print( LOG
"The number of characters differs within partition '$partition'! This partition has been omitted.\n"
                );
                return ();
            }
            $chars = $#aa;
            for ( my $j = $#aa ; $j >= 0 ; $j-- ) {
                $matrix->[$k][ $j + 1 ] = $aa[$j];
            }
        }
    }
    undef(@lines);

    ############################### MOP UNINFORMATIVE AND RECODE
    my @informative = ();
    my @numeric     = (
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B',
        'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
        'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'
    );
    for ( my $j = $#{ $matrix->[0] } ; $j > 0 ; $j-- ) {
        my %score  = ();
        my $states = 0;
        for ( my $k = $#{$matrix} ; $k >= 0 ; $k-- ) {
            if ( ( $matrix->[$k][$j] ne '?' ) && ( $matrix->[$k][$j] ne '-' ) )
            {
                if ( !exists( $score{ $matrix->[$k][$j] } ) ) {
                    $score{ $matrix->[$k][$j] } = 0;
                    $states++;
                }
                $score{ $matrix->[$k][$j] } += 1;
            }
        }
        if ( $states > 1 ) {
            my $minSteps = $states - 1;
            my $maxSteps = 0;
            my @steps    = sort( { $score{$b} <=> $score{$a} } keys(%score) );
            for ( my $i = $#steps ; $i > 0 ; $i-- ) {
                $maxSteps += $score{ $steps[$i] };
            }
            if ( ( $maxSteps - $minSteps ) > 0 ) {
                push( @informative, $j );
                if ($r) {
                    my %state = ();
                    $states = 0;
                    for ( my $k = $#{$matrix} ; $k >= 0 ; $k-- ) {
                        if (   ( $matrix->[$k][$j] ne '?' )
                            && ( $matrix->[$k][$j] ne '-' ) )
                        {
                            if ( $score{ $matrix->[$k][$j] } == 1 ) {
                                $matrix->[$k][$j] = '?';
                            }
                            else {
                                if ( !exists( $state{ $matrix->[$k][$j] } ) ) {
                                    $state{ $matrix->[$k][$j] } =
                                      $numeric[$states];
                                    $states++;
                                }
                                $matrix->[$k][$j] = $state{ $matrix->[$k][$j] };
                            }
                        }
                    }
                }
                if ( $maxStates < $states ) {
                    $maxStates = $states + 1;
                }
            }
        }
    }
    if ( $#informative >= 0 ) {
        $characters += $#informative + 1;
    }

    ############################### OUTPUT
    my $output = ();
    if ( $#informative >= 0 ) {
        if ($r) {
            $output = "&[num]\n";
        }
        else {
            $output = "&[prot nogaps]\n";
        }

        ############################### MAX LABEL SIZE
        my $max = 0;
        for ( my $k = $#{$matrix} ; $k >= 0 ; $k-- ) {
            if ( $max < length( $matrix->[$k][0] ) ) {
                $max = length( $matrix->[$k][0] );
            }
        }
        $max += 5;

        ############################### MATRIX
        for ( my $k = 0 ; $k <= $#{$matrix} ; $k++ ) {
            $output .= $matrix->[$k][0]
              . ( ' ' x ( $max - length( $matrix->[$k][0] ) ) );
            for ( my $j = $#informative ; $j >= 0 ; $j-- ) {
                $output .= $matrix->[$k][ $informative[$j] ];
            }
            $output .= "\n";
        }
        if ( $output ne "" ) { $pnames .= "$part: $partition\n" }

        #$output .= "\n";
    }
    $part++;
    if($output eq undef){print LOG "Partition $part: $partition is empty\n";}
	else{print LOG "Partition $part: $partition positions are [". join(",", @informative) ."]\n";}
    return ($output);
}
