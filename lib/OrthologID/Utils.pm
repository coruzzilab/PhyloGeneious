#
# OrthologID utilities
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
# Author: Ernest K Lee <elee@amnh.org>
#

package OrthologID::Utils;

use Text::Balanced qw( extract_bracketed );
use 5.008005;
use strict;
use warnings;

use Cwd;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
    'all' => [
        qw(

        )
    ]
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
  parseOIDgeneName
  readFasta
  writeFasta
  ungapAlignment
  untranslateTree
  getConsensus
  genSuperMatrix
  classifyQuery
  parentheticalTree
);

our $VERSION = '0.01';

#
# Given a #-delimited OrthologID gene name or fasta defline,
# return "species" and "gene/accession ID" portions
# as a list (species, geneID).  If there is no delimiter, undef is returned
# for the geneID.  Only the first # is a delimiter.  Everything after a
# whitespace is discarded.
#
sub parseOIDgeneName($) {
    my $defline = shift @_;
    my ( $species, $geneID );
    $defline =~ /^>?\s*([^#\s]*)#?(\S*)/;
    $species = $1;
    $geneID  = $2;
    return ( $species, $geneID );
}

#
# Return pointer to hash of sequences given list of FASTA file names.
# The filename "-" denotes standard input.
#
sub readFasta(@) {
    my @fastaFile = @_;
    my %seqs;
    foreach my $file (@fastaFile) {
        my $fh;
        if ( $file eq "-" ) {
            $fh = *STDIN;
        }
        else {
            open $fh, "$file"
              or warn "readFasta: Error opening FASTA file \"$file\": $!\n";
        }
        my ( $seqName, $seq ) = ( undef, undef );
        while (<$fh>) {
            chomp;
            next if /^\s*$/;       # skip blank lines
            if (/^>\s*(\S+)/) {    # defline
                                   # Save previous sequence
                $seqs{$seqName} = $seq if defined($seqName);
                $seq            = "";
                $seqName        = $1;
            }
            elsif ( defined $seq ) {
                s/[^A-Za-z\?\*-]+//g;    # remove unknown characters
                tr/a-z/A-Z/;             # store as uppercase
                $seq .= $_;
            }
            else {
                warn "readFasta: File not in FASTA format?\n";
                return undef;
            }
        }

        # Save last sequence
        $seqs{$seqName} = $seq if defined($seqName);

        close $fh if $file ne "-";
    }
    return \%seqs;
}

#
# Given pointer to sequence hash and file name, write sequences in fasta format
#
sub writeFasta($$) {
    my $seqRef  = shift;
    my $outFile = shift;

    open( FH, ">$outFile" ) or die "Cannot open $outFile for writing: $!\n";
    while ( my ( $name, $seq ) = each %$seqRef ) {
        print FH ">$name\n$seq";
    }
    close FH;
}

#
# Remove spanning gap regions from alignment.  These are columns of all gap or missing characters
# across all sequences.  This could happen when only a subset of sequences are chosen
# from a normal alignment.  Input is a pointer to sequence hash.
#
sub ungapAlignment($) {
    my $seqRef = shift;

    my @isGap;
    while ( my ( $name, $seq ) = each %$seqRef ) {

        #        if(defined $seq){print "$name was processed correctly\n";}
        #        else{print "$name has no sequence\n";}
        my @seqArr = split( //, $seq );
        foreach my $i ( 0 .. $#seqArr - 1 ) {
            if ( $seqArr[$i] ne '-' && $seqArr[$i] ne '?' ) {
                $isGap[$i] = 0;
            }
            elsif ( !defined( $isGap[$i] ) ) {
                $isGap[$i] = 1;
            }
        }
    }

    while ( my ( $name, $seq ) = each %$seqRef ) {
        my $newSeq = "";
        my @seqArr = split( //, $seq );
        foreach my $i ( 0 .. $#seqArr - 1 ) {
            $newSeq .= $seqArr[$i] if !$isGap[$i];
        }
        $seqRef->{$name} = $newSeq;
    }
}

#
# Given list of aligned sequences, return consensus sequence with differed positions
# replaced by ?'s.
#
sub getConsensus(@) {
    my @seqs = @_;

    my $len = length( $seqs[0] );
    foreach my $s (@seqs) {
        if ( length($s) != $len ) {
            warn "Sequences unaligned?";
            return undef;
        }
    }

    my @con = split( //, shift @seqs );    # consensus (as list of chars)

    foreach (@seqs) {
        my @list = split(//);
        foreach my $i ( 0 .. $#list ) {
            $con[$i] = '?' if $con[$i] ne $list[$i];
        }
    }
    return join( '', @con );
}

#
# Given a tree file, return list of untranslated parenthetical trees
#
sub untranslateTree($) {
    my $treeFile = shift;

    open FH, $treeFile or die "Cannot open tree file: $!";

    my $inTranslate = 0;    # inside translate block
    my @keyval;
    my @taxa;
    my @treeList;
    while (<FH>) {

        # Lines are not chomp'ed
        if (/^\s*Translate/i) {
            $inTranslate = 1;
        }
        elsif ($inTranslate) {
            if (/^\s*;/) { $inTranslate = 0; next }
            @keyval = split;
            chop( $keyval[1] ) if $keyval[1] =~ /[,;]$/;
            $taxa[ $keyval[0] ] = $keyval[1];
        }
        elsif (/^\s*tree[^\(]+(.*);/i) {
            ( my $untransTree = $1 ) =~ s/\d+/$taxa[$&]/ge;
            push @treeList, $untransTree;
        }
    }
    return @treeList;
}

sub parentheticalTree($) {
    my $treeFile = shift;
    my $tree;
    my @treeList;
    open FH, $treeFile or die "Cannot open tree file: $!";

    while (<FH>) {
        if (/^\(/) {
            $tree = $_;
            $tree =~ s/;//;
            $tree =~ s/\s//g;
            push @treeList, $tree;
        }
    }
    close FH;
    return @treeList;
}
#
# Change leading/trailing characters to the missing character.
# Input sequence string and missing char.  Return modified string.
#
sub gap2missing($$) {
    my $seq         = shift;
    my $missingChar = shift;
    my $tempStr;
    my $lenDiff;
    my $seqLen;

    $seqLen = length($seq);

    # replace leading gaps
    $tempStr = $seq;
    $tempStr =~ s/^-*//;
    $lenDiff = $seqLen - length($tempStr);
    $seq = $missingChar x $lenDiff . substr( $seq, $lenDiff ) if $lenDiff > 0;

    # replace trailing gaps
    $tempStr = $seq;
    $tempStr =~ s/-*$//;
    $lenDiff = $seqLen - length($tempStr);
    $seq     = substr( $seq, 0, length($tempStr) ) . $missingChar x $lenDiff
      if $lenDiff > 0;

    return $seq;

}

#
# Given lists of aligned sequences, returns a (concat) supermatrix as a string.
# Arguments: Pointer to list of taxon (species) names, pointer to list of
# outgroup taxa, pointer to list of (pipe-delimited)
# annotation for each partition (alignment), followed by list of pointers to
# hashes of alignments.  Missing taxa in each partition will be filled with
# unknown (?) characters.  Sequence names are assumed to start with a species
# prefix.  Only first of multiple sequences of same species will be used, no consensus.
#
sub genSuperMatrix {
    my $taxaRef     = shift;
    my $outgroupRef = shift;
    my $annoRef     = shift;
    my @taxa        = @$taxaRef;
    my %isOutgroup;
    my @sets        = ();
    my @setRange    = ();
    my $setNum      = 0;
    my $totLen      = 0;
    my $matrix      = "";
    my $missingChar = "?";    # change this if needed

    # init
    foreach (@$outgroupRef) {
        $isOutgroup{$_} = 1;
    }

    my $maxNameLen = 80;

    # Find max taxon name length for formatting
    foreach (@taxa) {
        my $len = length;
        $maxNameLen = $len if $len > $maxNameLen;
    }

    foreach (@_) {
        ungapAlignment($_);   # remove spanning gaps
        my %seqs = %$_;
        my $len = 0;
        my %curr;
        my $searchkey = (%seqs)[0];
        use List::MoreUtils 'first_index';
        my $ogindex = first_index { /$searchkey/i } @$annoRef;
        my @annokeys = split " ", (split /\|/, $$annoRef[$ogindex])[2];
        foreach my $name (@annokeys) {#while (my ($name, $seq) = each %seqs) {
            my $seq = $seqs{$name};
            
            # Calculate/Check length
            if ($len == 0) {
                $len = length($seq);
            }
            elsif ($len != length($seq)) {
                warn "Sequences not aligned!!\n";
                return undef;
            }
            
            # convert leading and trailing gaps to the missing char
            $seq = gap2missing($seq, $missingChar) if $seq !~ /^[^-?].*[^-?]$/;
            
            my $sp = (parseOIDgeneName($name))[0];
            if (defined($curr{$sp})) {
                #$curr{$sp} = getConsensus($curr{$sp}, $seq);
                next; # skip
            }
            else {
                $curr{$sp} = $seq;
            }
        }
        
        # Replace missing taxa with missing chars
        foreach my $taxon (@taxa) {
            if (! defined($curr{$taxon})) {
                $curr{$taxon} = $missingChar x $len;
            }
        }
        
        push(@sets, \%curr);
        
        # Record range
        $setRange[$setNum] = ($totLen+1) . "-" . ($totLen+$len);
        
        # Update running totals
        $totLen+= $len;
        $setNum++;
    }

    my @sizeSet;    # charsets by partition size (index), excluding outgroup

    $matrix .= "#NEXUS\n";
    $matrix .= "BEGIN DATA;\n";
    $matrix .= "DIMENSIONS NTAX=" . @taxa . " NCHAR=$totLen;\n";
    $matrix .=
"FORMAT INTERLEAVE SYMBOLS=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" DATATYPE=PROTEIN ";
    $matrix .= "MISSING=$missingChar GAP=-;\n";
    $matrix .= "MATRIX\n";
    foreach my $i ( 0 .. $setNum - 1 ) {
        my $anno = $$annoRef[$i];
        $anno =~ /^[^|]+/;
        my $partID =
          $&;  # assume first pipe-delimited field of annotation is partition id
               #my $partID = "OG" . ($i+1);
        $matrix .= "\n[ $anno|chars " . $setRange[$i] . " ]\n";
        my %currSet  = %{ $sets[$i] };
        my $partSize = 0;
        foreach my $taxon (@taxa) {
            $matrix .=
                "$taxon"
              . " " x ( $maxNameLen - length($taxon) + 1 )
              . $currSet{$taxon} . "\n";
            if ( $currSet{$taxon} !~ /^\?+$/ && !exists $isOutgroup{$taxon} ) {
                $partSize++;
            }

            #$partSize++ if $currSet{$taxon} !~ /^\?+$/;
        }

        # Add to charset for this size partition
        if ( defined $sizeSet[$partSize] ) {
            $sizeSet[$partSize] .= " $partID";
        }
        else {
            $sizeSet[$partSize] = " $partID";
        }

    }
    $matrix .= ";\nEND;\n\n";

    $matrix .= "BEGIN SETS;\n\n";
    foreach my $i ( 0 .. $setNum - 1 ) {
        my $partID = "OG" . ( $i + 1 );
        $matrix .= "CHARSET $partID = " . $setRange[$i] . ";\n";
    }
    $matrix .= "\n";

    foreach my $i ( 1 .. @taxa ) {
        $matrix .= "CHARSET SP$i =$sizeSet[$i];\n" if defined $sizeSet[$i];
    }
    $matrix .= "\n";

    $matrix .= "END;\n\n";

    return $matrix;
}

#
# Align and fit sequences to a profile (alignment)
# Arguments: sequences (to be align-fitted) and alignment files in FASTA format
# Return hash reference to profile-aligned (fitted) sequences
#
sub profileFit($$) {
    my $mafftProfileCmd = "mafft-profile";       # command must be in PATH
    my $seqFile         = shift;
    my $proFile         = shift;
    my $tmpFile         = "tmp.fa";
    my $noMissingFile   = "profile-nomiss.fa";
    my %pafSeq;

    my $seqRef = readFasta($seqFile);

# Make sure profile has no missing (?) characters -- mafft-profile ignores them.
# Otherwise, replace with '-'.
    open( PFH, $proFile ) or die "Cannot open profile $proFile: $!";
    my $savRS = $/;
    undef $/;
    my $profile   = <PFH>;
    my $noMissing = $profile;
    $noMissing =~ tr/?/-/;
    if ( $noMissing ne $profile ) {
        open( NFH, ">$noMissingFile" )
          or die "Cannot open $noMissingFile for writing: $!";
        print NFH "$noMissing";
        close NFH;
        $proFile = $noMissingFile;
    }
    $/ = $savRS;
    close PFH;

    # Write out each sequence and profile align-fit
    while ( my ( $name, $seq ) = each %$seqRef ) {
        writeFasta( { $name => $seq }, $tmpFile );

        # Read output alignment into single string
        $savRS = $/;
        undef $/;
        my $output = `$mafftProfileCmd $tmpFile $proFile 2>/dev/null`;
        $/ = $savRS;

        # Store alignment sequences into array
        $output =~ s/>.*/>/g;
        $output =~ s/\s//g;
        my @seqs = split( />/, $output );

        # Transpose
        my @cPos;    # character position
        foreach my $n ( 0 .. $#seqs ) {
            my @chars = split( //, $seqs[$n] );
            foreach my $k ( 0 .. $#chars ) {
                $cPos[$k] .= $chars[$k];
            }
        }

        # Remove characters that cause gaps in original matrix
        my $fittedSeq = "";
        foreach (@cPos) {
            $fittedSeq .= substr( $_, 0, 1 ) if !/^.\-+$/;
        }

        $pafSeq{$name} = $fittedSeq;
    }
    return \%pafSeq;
}

#
# Classify query sequences into tree, given matrix (alignment and tree).
# Arguments: alignment FASTA file, (unaligned) query FASTA file, guide tree file, and
#            reference to list of outgroup
# Return parenthetical tree with query in it.
# If the outgroup list is empty, the tree is unrooted.
# If the guide tree is uninformative, thus cannot be used as a backbone constraint,
# swapping will be done.
#
# *** OBSOLETE ***
#
sub classifyQueryBackboneConstraint($$$$) {
    my $alignFile     = shift;
    my $queryFile     = shift;
    my $guideTreeFile = shift;
    my @outgroup      = @{ $_[0] };
    my $treeFile      = "classify.tre";
    my $uninf;    # if guide tree is informative

    my $paupCmd = "paup";    # must be in PATH

    # Check number of taxa
    my @taxa = ();
    open FH, $alignFile or die "Cannot open alignment file: $!";
    while (<FH>) {
        chomp;
        push @taxa, $1 if /^>(\S+)/;
    }
    close FH;

#    if (@taxa < 4) {
#        $uninf = 1;
#        print STDERR "guide tree with less than 4 taxa ... will not use contraint\n";
#    }
    open FH, $queryFile or die "Cannot open query file: $!";
    while (<FH>) {
        chomp;
        push @taxa, $1 if /^>(\S+)/;
    }
    close FH;
    if ( @taxa < 4 ) {    # nothing to do
        my $t = join ',', @taxa;
        print STDERR "  ==> ($t) returned\n";
        return "($t)";
    }

    #    if (!defined($uninf)) {
    # Check guide tree
    open GFH, $guideTreeFile or die "Cannot open guide tree file: $!";
    while (<GFH>) {
        chomp;
        next if !/^\s*Tree\s*\S+\s*=[^(]*(.*);/i;    # tree line
        my $parenTree = $1;
        if ( $parenTree =~ /^(\(\([^)]+\),[^,)]+\)|\([^,)]+,\([^)]+\)\))$/ ) {
            print STDERR
"guide tree $parenTree equivalent to all trees ... will not use constraint\n";
            $uninf = 1;
        }
        s/[^(]//g;
        if ( length == 1 ) {
            print STDERR
"guide tree $parenTree uninformative ... will not use contraint\n";
            $uninf = 1;
        }
    }

    #    }

    my $qSeqRef  = profileFit( $queryFile, $alignFile );
    my $seqRef   = readFasta($alignFile);
    my @qNames   = keys %$qSeqRef;
    my @seqNames = keys %$seqRef;
    my $nSeq     = @qNames + @seqNames;
    my $len      = length $qSeqRef->{ $qNames[0] };

    # Create matrix
    my $nexFile = "classify.nex";
    open( NEX, ">$nexFile" ) or die "Failed to open $nexFile: $!\n";
    print NEX "#NEXUS\nbegin data;\n";
    print NEX "dimensions ntax=$nSeq nchar=$len;\n";
    print NEX "format missing=? symbols=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" ";
    print NEX "datatype=PROTEIN gap=-;\n\n";
    print NEX "matrix\n";

    while ( my ( $taxon, $seq ) = each %$seqRef ) {
        print NEX "$taxon\t$seq\n";
    }
    foreach (@qNames) {
        print NEX "$_\t$qSeqRef->{$_}\n";
    }
    print NEX "\n;\nend;\n";
    print NEX "begin paup;\n";
    print NEX "set increase=auto notifybeep=no;\n";
    print NEX "outgroup @outgroup;\n" if @outgroup;
    if ( defined $uninf ) {    # uninformative guide tree, do swap
        if ( @taxa <= 10 ) {    # exact search
            print NEX "bandb;\n";
        }
        else {
            print NEX "hsearch addseq=random swap=tbr timelimit=300;\n";
        }
    }
    else {    # use guide tree as backbone
        print NEX "loadconstr file=$guideTreeFile asbackbone;\n";
        if ( @taxa <= 10 ) {    # exact search
            print NEX "bandb enforce=yes;\n";
        }
        else {
            # Also do NNI swapping if more than 1 query
            if ( @qNames > 1 ) {
                print NEX
"hsearch addseq=random enforce=yes swap=nni rearrlimit=10000 limitperrep=yes;\n";
                print NEX "filter best;\n";
                print NEX "contree all/file=tmp.tre replace=yes;\n";
            }
            print NEX "hsearch addseq=closest enforce=yes swap=none;\n";
            if ( @qNames > 1 ) {
                print NEX "gettrees file=tmp.tre warntree=no mode=7;\n";
                print NEX "filter best;\n";
            }
        }
    }

    # Save trees and score of first tree
    print NEX "pscores /scorefile=tree-score replace=yes;\n";
    print NEX "savetrees file=$treeFile replace=yes";
    print NEX " root=yes" if @outgroup;
    print NEX ";\nquit;\nend;\n";

    # Run PAUP*
    my $status = system("$paupCmd -n $nexFile >/dev/null");
    if ( $status != 0 ) {
        die "PAUP did not execute successfully: $?";
        return undef;
    }

    my @tree = untranslateTree($treeFile);
    return $tree[0];    # return first tree
}

#
# Generate all possible trees given constraint tree by adding taxon
# to all possible positions.
#
# Params: parenthetical tree, new taxon name
#
sub genAllTrees($$) {
    my $parenTree = shift;
    my $newTaxon  = shift;
    my @allTrees  = ();

    my @chars = split( / */, $parenTree );

    # Pair with all taxa
    while ( $parenTree =~ /[\w\d#]+/g ) {
        push( @allTrees, "$`($&,$newTaxon)$'" );
    }

    # Add to all clades (polytomies)
    while ( $parenTree =~ /\)/g ) {
        push( @allTrees, "$`,$newTaxon$&$'" );
    }

    # Pair with all clades
    while ( $parenTree =~ /\(/g ) {
        my ( $head,    $tail )   = ( $`, "$&$'" );
        my ( $matched, $suffix ) = extract_bracketed( $tail, '()' );
        push( @allTrees, "$head($newTaxon,$matched)$suffix" );
    }

    return @allTrees;
}

#
# Classify query sequences into tree, given matrix (alignment and tree), by placing
# sequences into all possible positions of the tree, and return the shortest tree.
#
# Arguments: alignment FASTA file, (unaligned) query FASTA file, constraint/guide
#            tree file, and reference to list of outgroup
#
# Return value: parenthetical tree with query in it
#
# If the outgroup list is empty, the tree is unrooted.
#
# ** need to test multi-fasta query file
#
#
sub classifyQuery($$$$) {
    my $alignFile     = shift;
    my $queryFile     = shift;
    my $guideTreeFile = shift;
    my @outgroup      = @{ $_[0] };
    my $bestTreeFile  = "classify.tre";

    my $paupCmd = "paup";    # must be in PATH

    # Get constraint tree taxa
    my @taxa = ();
    my %taxaSeen;
    open FH, $alignFile or die "Cannot open alignment file: $!";
    while (<FH>) {
        chomp;
        if (/^>(\S+)/) {

            # Normally, taxon names should not be repeated.
            # This guards against duplicates in OID input data.
            if ( !defined( $taxaSeen{$1} ) ) {
                push @taxa, $1;
                $taxaSeen{$1} = 1;
            }
        }
    }
    close FH;

    # Get queries
    my @qTaxa = ();
    open FH, $queryFile or die "Cannot open query file: $!";
    while (<FH>) {
        chomp;
        push @qTaxa, $1 if /^>(\S+)/;
    }
    close FH;

    my $totTaxa            = @taxa;
    my $unalignedQuery     = readFasta($queryFile);
    my $augmentedAlignFile = $alignFile;
    my $tempTreeFile       = $guideTreeFile;
    foreach my $q ( keys(%$unalignedQuery) ) {

        $totTaxa++;

        # Profile align query into matrix and get new number of chars
        my $singleQueryFile = "single_query";
        open FH, ">$singleQueryFile"
          or die "Cannot open \"$singleQueryFile\" for writing: $!\n";
        print FH ">$q\n";
        print FH $unalignedQuery->{$q} . "\n";
        close FH;
        my $qSeqRef = profileFit( $singleQueryFile, $augmentedAlignFile );
        my $seqRef  = readFasta($augmentedAlignFile);

        $augmentedAlignFile =
          "tmp.aln";    # alignment file augmented with this query
        open AFH, ">$augmentedAlignFile"
          or die "Cannot open \"$augmentedAlignFile\" for writing: $!\n";
        while ( my ( $r, $s ) = each %$seqRef ) {
            print AFH ">$r\n$s\n";
        }
        print AFH ">$q\n$qSeqRef->{$q}\n";
        close AFH;

        my @seqNames = keys %$seqRef;
        my $len      = length $qSeqRef->{$q};

        #unlink $singleQueryFile;

        # Create matrix
        my $nexFile = "classify.nex";
        open( NEX, ">$nexFile" ) or die "Failed to open $nexFile: $!\n";
        print NEX "#NEXUS\nbegin data;\n";
        print NEX "dimensions ntax=$totTaxa nchar=$len;\n";
        print NEX "format missing=? symbols=\"ABCDEFGHIKLMNPQRSTUVWXYZ\" ";
        print NEX "datatype=PROTEIN gap=-;\n\n";
        print NEX "matrix\n";

        while ( my ( $taxon, $seq ) = each %$seqRef ) {
            print NEX "$taxon\t$seq\n";
        }
        print NEX "$q\t$qSeqRef->{$q}\n";
        print NEX "\n;\nend;\n";
        print NEX "begin paup;\n";
        print NEX "set increase=auto;\n";
        print NEX "end;\n";

        # Generate all possible trees
        my $translateBlk = "";
        my $rootFlag;
        my $parenTree;
        print NEX "begin trees;\n";
        open GFH, $tempTreeFile
          or die "Cannot open guide tree \"$tempTreeFile\": $!\n";
        while (<GFH>) {
            if (/^\s*Translate/i) {
                $translateBlk .= "$_\t\t$totTaxa $q,\n";
            }
            elsif (/^\s*\d+\s+\w+/) {
                $translateBlk .= $_;
            }
            elsif (/^tree PAUP_1 = (\[\&(?:U|R)\]) (\(.*\));/i) {
                $rootFlag  = $1;
                $parenTree = $2;
            }
        }
        print NEX "$translateBlk;\n";
        close GFH;

        my @allTrees = genAllTrees( $parenTree, $totTaxa );
        my $cnt      = 1;
        foreach my $t (@allTrees) {
            print NEX "tree CLASSIFY_" . $cnt++ . " = $rootFlag $t;\n";
        }
        print NEX "end;\n";

        # Find shortest tree
        my $scoreFile = "tree-score";
        my ( $shortestTreeNum, $shortestTree );
        print NEX "begin paup;\n";
        print NEX "filter best;\n";
        print NEX "pscores /scorefile=$scoreFile replace=yes;\n";
        print NEX "quit;\nend;\n";
        close NEX;

        my $status = system("$paupCmd -n $nexFile >/dev/null");
        if ( $status != 0 ) {
            die "PAUP did not execute successfully: $?";
            return undef;
        }

        open SFH, $scoreFile or die "Cannot open score file: $!\n";
        <SFH>;
        my $line = <SFH>;
        chomp $line;
        $shortestTreeNum = ( split( /\s+/, $line ) )[0];
        $shortestTree    = $allTrees[ $shortestTreeNum - 1 ];
        close SFH;

        # Write out the (translated) tree
        open TRE, ">$bestTreeFile" or die "Cannot open \"$bestTreeFile\": $!\n";
        print TRE "#NEXUS\n";
        print TRE "Begin trees;\n";
        print TRE "$translateBlk;\n";
        print TRE "tree CLASSIFY_$shortestTreeNum = $rootFlag $shortestTree;\n";
        print TRE "End;\n";
        close TRE;

        $tempTreeFile = $bestTreeFile;
    }

    my @tree = untranslateTree($bestTreeFile);
    return $tree[0];    # return tree
}
