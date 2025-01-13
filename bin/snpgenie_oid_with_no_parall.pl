#! /usr/bin/perl
## /usr/local/software/PERL/perl-5.26.0/bin/perl

# PROGRAM: SNPGenie Nei-Gojobori algorithm adapted for OrthologID and HPC

# Copyright (C) 2017 Chase W. Nelson

#########################################################################################
### NOTES TO KRANTHI & GIL:
### (1) a head script should call this once for a single $OID_USER_DIR/data/<FAMILY> directory.
### (2) results are placed in files in that family directory.
### (3) perhaps use partitionMembers file to determine which families to call, since not
###     all families contain partitions. Results files will be empty if no partitions.
### (4) all partitions in the family will be sequentially processed. Thus, number of partitions
###     in family is proportional to runtime.
### (5) some parallelism will be invoked using the Perl module Parallel::ForkManager, controlled
###     by the --procs_per_node argument. This will parallelize the treatment of codons at
###     certain points in the code, but NOT the processing of partitions.
#########################################################################################

#########################################################################################
### EXAMPLE CALL:
# $OID_HOME/bin/snpgenie_oid_HPC.pl --family_directory=<dir_name> --ignore_done_flag --between_group_flag --procs_per_node=<integer> --sliding_window_size=<integer_gt_0> --sliding_window_step=<integer_gt_0> --num_bootstraps=<integer_gt_1>
# $OID_HOME/bin/snpgenie_oid_HPC.pl --family_directory=10 --ignore_done_flag --between_group_flag --procs_per_node=1 --sliding_window_size=10 --sliding_window_step=1 --num_bootstraps=1000
#########################################################################################

#########################################################################################
## LICENSE
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################################

# DATE CREATED: January 18, 2017

# AUTHOR: Chase W. Nelson

# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com

# AFFILIATION1: Sackler Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA
# AFFILIATION2: Special Volunteer, Division of Cancer Epidemiology & Genetics, National
#     Cancer Institute, National Institutes of Health, Rockville, MD 20850, USA
# AFFILIATION3: BigPlant Consortium, Center for Genomics and Systems Biology, New York
#     University, New York, NY 10003, USA

# CITATION1: SNPGenie, https://github.com/chasewnelson/snpgenie
# CITATION2: Nelson CW, Moncla LH, Hughes AL (2015) SNPGenie: estimating evolutionary
#	parameters to detect natural selection using pooled next-generation sequencing data.
#	Bioinformatics 31(22):3709-11, doi: 10.1093/bioinformatics/btv449.

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;

#use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw(max);

# Make sure we have (and store) the environment variables
my ( $OID_HOME, $OID_USER_DIR );

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
}

STDOUT->autoflush(1);

#########################################################################################
# INITIALIZE VARIABLES (all optional and/or defaulted)
my $family_directory;
my $ignore_done_flag; # ignore .done, auto-remove all prior results files, re-do
my $between_group_flag;

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions(
    "ignore_done_flag" =>
      \$ignore_done_flag,    # optional Boolean; set to false if not given
    "between_group_flag" =>
      \$between_group_flag,    # optional Boolean; set to false if not given
    "family_directory=s" => \$family_directory    # mandatory string parameter
  )

  or die
"\n### WARNING: Error in command line arguments. SNPGenie for OID terminated.\n\n";

# N.B.: When an argument, e.g., sliding_window_size, is called only as a flag, its value is 0
# When it is not called at all, it is null

# Get the time
my $time1       = time;
my $local_time1 = localtime;

print
"\n#########################################################################################\n";
print "SNPGenie for Phylogeneous initiated at local time $local_time1\n";

# Test and/or reset OPTIONS given the user's INPUT and RECORD PARAMETERS
if ( !$family_directory ) {
    die "\n### WARNING: The --family_directory option must be specified\n"
      . "### SNPGenie for OID terminated.\n\n";
}

if ( !$ignore_done_flag ) {
    $ignore_done_flag = 0
      ; # default behavior: re-do nothing if flag absent, i.e., pick up where left off
}
else {
    $ignore_done_flag = 1;
}

if ( !$between_group_flag ) {
    $between_group_flag = 0;    # default behavior: no between-group analysis
}
else {
    $between_group_flag = 1;
}

my $family_id = $family_directory;

if ( $family_directory =~ /(\d+)/ ) {
    $family_id = $1;
}

print
"#########################################################################################\n";
print "\nWorking on family $family_id\n";

if ( $ignore_done_flag == 1 ) {
    unlink "$OID_USER_DIR/data/$family_id/within_group_results.txt";
    unlink "$OID_USER_DIR/data/$family_id/within_group_codon_results.txt";
    unlink "$OID_USER_DIR/data/$family_id/this_family_log.txt";
    unlink "within\_group\_codon\_results\.txt";
    unlink "$OID_USER_DIR/data/$family_id/retained_aa_positions.txt";

    #	unlink "$OID_USER_DIR/data/$family_id/aligned.cds";
    unlink "$OID_USER_DIR/data/$family_id/.snpgenie_done";
    unlink "$OID_USER_DIR/data/$family_id/.snpgenie_ignore";

#	unlink "$OID_USER_DIR/data/$family_id/family$family_id\_between_group_codon_results.txt";
#	unlink "$OID_USER_DIR/data/$family_id/family$family_id\_between_group_results.txt";
    unlink "$OID_USER_DIR/data/$family_id/between_group_codon_results.txt";
    unlink "$OID_USER_DIR/data/$family_id/between_group_results.txt";

  #	unlink "between\_group\_sw_" . $sliding_window_size . "codons\_results.txt";
###	unlink "$OID_USER_DIR/data/$family_id/family$family_id\_between_group_polymorphism.txt";
}

if ( -f "$OID_USER_DIR/data/$family_id/.snpgenie_done" )
{    # if we're ignoring done files, we just deleted this
    print "Already completed. Moving on to next family.\n";

#next FAMILY; # this doesn't work because a child process apparently cannot be called
}
else {

    # There's not a done file, so let's just delete if they're there to be sure
    unlink "$OID_USER_DIR/data/$family_id/within_group_results.txt";
    unlink "$OID_USER_DIR/data/$family_id/within_group_codon_results.txt";
    unlink "$OID_USER_DIR/data/$family_id/this_family_log.txt";
    unlink "within\_group\_codon\_results\.txt";
    unlink "$OID_USER_DIR/data/$family_id/retained_aa_positions.txt";

    #	unlink "$OID_USER_DIR/data/$family_id/aligned.cds";
    unlink "$OID_USER_DIR/data/$family_id/.snpgenie_done";
    unlink "$OID_USER_DIR/data/$family_id/.snpgenie_ignore";
    unlink "$OID_USER_DIR/data/$family_id/between_group_codon_results.txt";
    unlink "$OID_USER_DIR/data/$family_id/between_group_results.txt";

  #	unlink "between\_group\_sw_" . $sliding_window_size . "codons\_results.txt";
###	unlink "$OID_USER_DIR/data/$family_id/family$family_id\_between_group_polymorphism.txt";

    # Go to the family directory
    chdir("$OID_USER_DIR\/data\/$family_id");

    my @matching_files = glob "*.aligned";

    # There is a .aligned file for the FAMILY
    if ( scalar(@matching_files) == 1 ) {

        #print "\nWe've got one matching file\n";
        my $this_family_file_name = $matching_files[0];

        ##################################################################################
        # STORE INCLUDED PARTITIONS and determine WHICH CDS are needed
        my %family_partition_members_hh;
        my %family_members_h;
        my %family_species_h;
        my %family_partition_species_hh;

        # For example, keys of family_members_h might be:
        #AaphrophilusNJ8700#YP_003007290.1
        #Msuccinici#YP_089554.1
        #Asuccinogenes130Z#YP_001344574.1
        #HparainfATCC33392#HMPREF9417_0308
        #EcoliK12MG1655#NP_417469.1
        #Apleuropneumoniae10D13039#appser10_13970

        chdir("$OID_USER_DIR");

        if ( -f "$OID_USER_DIR/partitionMembers.3.txt" ) {
            open( IN_PARTITIONS, "partitionMembers.3.txt" )
              or die "Could not open file partitionMembers.3.txt\n";

            print
"Recording family/partition members from partitionMembers.3.txt...\n";

# Example file lines, partitionMembers.txt
#p1 = Arath#AT1G05700.1 Glyma#Glyma07g01620.1 	Family1
#p2 = Arath#AT1G51620.2 Glyma#Glyma01g00790.1 Teich#Gene.11002 Grobu#g.89669 Solyc#Solyc03g093380.1.1 	Family1

# EXAMPLE & EXPLANATION:
# Arath#AT1G05700.1 can be found in $OID_USER_DIR/blastdb within files Arath_cds.fa (DNA) and/or Arath.faa (AA)
# within the cds file, we have the header >Arath#AT1G05700.1
# these headers will be stored in @members_array

            while (<IN_PARTITIONS>) {
                chomp;

               # if(/^(\d+?): \[ \w+\|\w+\|([\w\#\.\s]+)\|chars \d+\-\d+ \]$/) {
                if (/^p(\d+) = ([\w\#\.\s]+)\s+Family(\d+)$/) {

                    my $partition_number = $1;

                    # print "\nPartition number: $partition_number\n";

                    my $members_string = $2;

                    # print "\nMembers string: $members_string\n";

                    my $family_number = $3;

                    if ( $family_number == $family_id ) {

                        # print "\nFamily number: $family_number\n";

                        my @members_array = split( /\s+/, $members_string );

                        # print "\nMembers array: @members_array\n";
                        foreach my $this_member (@members_array) {
                            $family_partition_members_hh{$family_number}
                              ->{$partition_number}->{$this_member} = 1;
                            print
"Added $this_member to family $family_number partition $partition_number\n";
                            $family_members_h{$this_member} = 1;

                            # print "\n\nSTORING member $this_member\n";

                            if ( $this_member =~ /(\w+)#/ ) {

                                # print "\n\nSTORING species $1\n";
                                my $this_species = $1;

                                $family_species_h{$this_species} = 1;

                                if (
                                    $family_partition_species_hh{
                                        $partition_number}->{$this_species} == 1
                                  )
                                {
                                    warn
"### WARNING: there are multiple members of species $this_species "
                                      . "in family $family_number partition $partition_number.\n"
                                      . "### Please fix your data. There should be one representative ortholog for each species in a partition.\n"
                                      . "### One gene will be arbitrary chosen to proceed.\n\n";
                                }

                                $family_partition_species_hh{$partition_number}
                                  ->{$this_species} = 1;
                            }

                        }
                    }
                }
            }

            close IN_PARTITIONS;
        }
        elsif ( -f "$OID_USER_DIR/partitionMembers.txt" ) {
            open( IN_PARTITIONS, "partitionMembers.txt" )
              or die "Could not open file partitionMembers.txt\n";

            print
"Recording family/partition members from partitionMembers.txt...\n";

# Example file lines, partitionMembers.txt
#p1 = Arath#AT1G05700.1 Glyma#Glyma07g01620.1 	Family1
#p2 = Arath#AT1G51620.2 Glyma#Glyma01g00790.1 Teich#Gene.11002 Grobu#g.89669 Solyc#Solyc03g093380.1.1 	Family1

# EXAMPLE & EXPLANATION:
# Arath#AT1G05700.1 can be found in $OID_USER_DIR/blastdb within files Arath_cds.fa (DNA) and/or Arath.faa (AA)
# within the cds file, we have the header >Arath#AT1G05700.1
# these headers will be stored in @members_array

            while (<IN_PARTITIONS>) {
                chomp;

               # if(/^(\d+?): \[ \w+\|\w+\|([\w\#\.\s]+)\|chars \d+\-\d+ \]$/) {
                if (/^p(\d+) = ([\w\#\.\s]+)\s+Family(\d+)$/) {

                    my $partition_number = $1;

                    # print "\nPartition number: $partition_number\n";

                    my $members_string = $2;

                    # print "\nMembers string: $members_string\n";

                    my $family_number = $3;

                    if ( $family_number == $family_id ) {

                        # print "\nFamily number: $family_number\n";

                        my @members_array = split( /\s+/, $members_string );

                        # print "\nMembers array: @members_array\n";
                        foreach my $this_member (@members_array) {
                            $family_partition_members_hh{$family_number}
                              ->{$partition_number}->{$this_member} = 1;
                            print
"Added $this_member to family $family_number partition $partition_number\n";
                            $family_members_h{$this_member} = 1;

                            # print "\n\nSTORING member $this_member\n";

                            if ( $this_member =~ /(\w+)#/ ) {

                                # print "\n\nSTORING species $1\n";
                                my $this_species = $1;

                                $family_species_h{$this_species} = 1;

                                if (
                                    $family_partition_species_hh{
                                        $partition_number}->{$this_species} == 1
                                  )
                                {
                                    warn
"### WARNING: there are multiple members of species $this_species "
                                      . "in family $family_number partition $partition_number.\n"
                                      . "### Please fix your data. There should be one representative ortholog for each species in a partition.\n"
                                      . "### One gene will be arbitrary chosen to proceed.\n\n";
                                }

                                $family_partition_species_hh{$partition_number}
                                  ->{$this_species} = 1;
                            }

                        }
                    }
                }
            }

            close IN_PARTITIONS;
        }
        elsif ( -f "$OID_USER_DIR/includedPartitions\.3\.txt" ) {

            my $included_partitions_file =
              "$OID_USER_DIR/includedPartitions\.3\.txt";

            open( IN_PARTITIONS, "$included_partitions_file" )
              or die "Could not open file $included_partitions_file\n";

            print
"Recording family/partition members from $included_partitions_file...\n";

# Example file lines, included.partitions.2.txt
#2: [ OG2|Family1|Arath#AT1G51620.2 Glyma#Glyma01g00790.1 Glyma#Glyma07g15270.1 Glyma#Glyma07g15270.2 Teich#Gene.11002 Grobu#g.89669 Grobu#g.89670 Solyc#Solyc03g093380.1.1 Solyc#Solyc03g121230.2.1|chars 1346-2439 ]
#3: [ OG3|Family1|Teich#Gene.10387 Arath#AT1G67720.1 Teich#Gene.6530 Betav#Bv6_146480_fkrc.t1 Grobu#g.38678 Solyc#Solyc05g014240.2.1 Glyma#Glyma08g10640.1 Glyma#Glyma05g27650.1 Glyma#Glyma05g27650.2 Glyma#Glyma18g01450.1 Glyma#Glyma11g37500.3 Glyma#Glyma11g37500.1 Glyma#Glyma11g37500.2|chars 2440-3405 ]

# Example file lines, included.partitions.3.txt
#1	[ OG1|Family10|Pquad#c9235 Btola#c20757 Btola#c20755 Btola#c20756 Bboli#c19211 Ccana#c21197 Scana#c21871 Scana#c19237 Grobu#c24617 Ftrin#c8635 Tmult#c10387 Eangu#c23591 Atrif#c18316 Atrif#c18317 Aarte#c13199 Aarte#c13198 Cpala#c6472 Spuch#c32360 Schry#c8702|chars 1-497 ]
#2	[ OG2|Family10|Arath#AT4G24400.1 Teich#c10297 Fchil#c19313 Fchil#c9559 Fchil#c9560 Clute#c16740 Glyma#Glyma02g38180.2 Clute#c16739 Among#c19871 Langu#c24825 Loreo#c11264 Lsubi#c22894 Lsubi#c22893 Lsubi#c22895 Lsubi#c22896 Ahypo#c15195 Aspin#c20109 Aerin#c24200 Aspin#c9261 Glyma#Glyma19g05410.3 Glyma#Glyma06g09700.2 Glyma#Glyma04g09610.1 Bdomi#c18556 Ainca#c17840 Ainca#c17839 Svulg#c19062 Pbryo#c21963 Aimbr#c13165 Arose#c10979 Cconm#c11334 Cconm#c11335 Pminu#c11568 Cconm#c24009 Aadsc#c20295 Catac#c12252 Clana#c7663 Cpala#c12846 Eangu#c13090 Atrif#c9788 Aarte#c13531 Ftrin#c18674 Tmult#c22285 Eangu#c18154 Grobu#c33565 Bboli#c23282 Btola#c19445 Pquad#c20778 Ccana#c21216 Scana#c23012 Schry#c10112 Spuch#c29056 Schry#c7887 Schry#c15758 Aspat#c18058 Mmono#c18111 Mcras#c7093 Casia#c13545 Adese#c6964 Pdulc#c7638 Sinca#c13441 Sinca#c16794 Einte#c20999 Einte#c21000 Einte#c20998 Dmete#c7883 Fdenu#c19896 Pinte#c15082 Dmete#c7882 Schil#c18330 Solyc#Solyc04g076810.2.1 Ppinn#c21127 Pcamp#c17306|chars 498-1033 ]
#3	[ OG3|Family10|Zjapo#c10011 Zjapo#c10010 Aadsc#c15570 Sgran#c8970 Nnard#c16801 Nnard#c16802 Jfrig#c26807 Sgran#c7225 Bradi#Bradi1g36520.1 Bsylv#c22584 Bsylv#c22583 Dangu#c9234 Fprat#c21569 Fchry#c20122 Ccurv#c23417 Ccurv#c23416 Eglab#c4780 Zeama#GRMZM2G137569_T01 Zeama#GRMZM2G171507_T02 Bdact#c25951 Bsimp#c18286 Mdecu#c17267 Mdecu#c17268 Mdecu#c17266|chars 1034-1561 ]

            while (<IN_PARTITIONS>) {
                chomp;
                if (
/^(\d+?)\s+\[\s+\w+\|Family(\d+)\|([\w\#\.\s]+)\|chars \d+\-\d+\s+\]$/
                  )
                {

                    my $partition_number = $1;

                    print "\nPartition number: $partition_number\n";

                    my $members_string = $3;

                    print "\nMembers string: $members_string\n";

                    my $family_number = $2;

                    if ( $family_number == $family_id ) {
                        print "\nFamily number: $family_number\n";

                        my @members_array = split( /\s+/, $members_string );

                        print "\nMembers array: @members_array\n";
                        foreach my $this_member (@members_array) {
                            $family_partition_members_hh{$family_number}
                              ->{$partition_number}->{$this_member} = 1;
                            print
"Added $this_member to family $family_number partition $partition_number\n";

                            $family_members_h{$this_member} = 1;
                            print "\n\nSTORING member $this_member\n";

                            if ( $this_member =~ /(\w+)#/ ) {
                                print "\n\nSTORING species $1\n";
                                my $this_species = $1;

                                $family_species_h{$this_species} = 1;

                                if (
                                    $family_partition_species_hh{
                                        $partition_number}->{$this_species} == 1
                                  )
                                {
                                    warn
"### WARNING: there are multiple members of species $this_species "
                                      . "in family $family_number partition $partition_number.\n"
                                      . "### Please fix your data. There should be one representative ortholog for each species in a partition.\n"
                                      . "### One gene will be arbitrary chosen to proceed.\n\n";
                                }

                                $family_partition_species_hh{$partition_number}
                                  ->{$this_species} = 1;
                            }

                        }
                    }
                }

            }

            close IN_PARTITIONS;

        }    # done with partitionMembers

        #		print "\n\nKeys of family_members_h:\n";
        #		foreach(keys %family_members_h) {
        #			print "$_\n";
        #		}

        # Get the .cds sequence file names and stay in the /blastdb/ directory
        chdir("$OID_USER_DIR\/blastdb");

        #my @cds_file_names = glob "*.cds";
        my @cds_file_names = glob "*_cds.fa";
        my @cds_to_add     = glob "*.cds";
        @cds_file_names = ( @cds_file_names, @cds_to_add );

#unless(scalar(@OID_USER_DIR_contents) > 0) { die "\n\n# There are no directories in $OID_USER_DIR\/data\n\n"; }

        unless ( scalar(@cds_file_names) > 0 ) {
            die
"\n\n# There are no .cds nucleotide files in $OID_USER_DIR\/blastdb\n\n";
        }

# READ IN ONLY those .cds sequences from the NUCLEOTIDE fasta files that are NEEDED
        my %header_to_index_hash
          ;    # this may later be problematic because of memory usage
        my $seq = '';
        my @seqs_arr;
        my $header = '';
        my @headers_arr;
        my $seq_num = 0;

        print "\n";

        my $add_current_flag = 0;

        ##################################################################################
        # Store needed CDS data
        my $seen_match = 0;

        foreach my $cds_file (@cds_file_names) {

            open( IN_FASTA, "$cds_file" )
              or die "Could not open file $cds_file\n";

            # print "Recording coding sequence data for file $cds_file...\n";

            while (<IN_FASTA>) {
                chomp;

                my $this_line = $_;

                if ( $this_line =~ />/ ) {
                    my $this_line_no_carrot = $this_line;
                    $this_line_no_carrot =~ s/>//;

# print "\nthis_line is: $this_line\nthis_line_no_carrot is: $this_line_no_carrot\n";

                    if ( $family_members_h{$this_line_no_carrot} == 1 ) {

                # print "We've matched family member $this_line_no_carrot...\n";
                        $add_current_flag = 1;
                        $seen_match       = 1;

                        print
"Recording nucleotide sequence $this_line_no_carrot from file $cds_file...\n";

                        if ( $seq_num == 0 ) {

                            #$header = $_;
                            $header = $this_line;
                            $seq_num++;
                        }
                        else {
                            $seq = uc($seq);
                            $seq =~ tr/U/T/;
                            push( @seqs_arr,    $seq );
                            push( @headers_arr, $header );

           # print "\nStoring in index ". ($seq_num - 1) ." $header...\n$seq\n";
           # print "\nseq_num $seq_num header $header: $seq\n";

                            $header_to_index_hash{$header} = ( $seq_num - 1 );

                            #$header = $_;
                            $header = $this_line;
                            $seq_num++;

                            $seq = '';
                        }
                    }
                    else
                    {   # else don't store current, since it's not in the family

                        if ( $add_current_flag == 1 ) {
                            $add_current_flag = 0;
                        }
                    }
                }
                else {    # not a header line
                    if ( $add_current_flag == 1 ) {

                        # $seq .= $_;
                        $seq .= $this_line;
                    }     # else not storing, since it's not in the family
                }
            }

            close IN_FASTA;
        }    # last CDS file

        # Store last one
        if ( $seen_match == 1 ) {
            if ( !$header_to_index_hash{$header} ) {
                $seq = uc($seq);
                $seq =~ tr/U/T/;
                push( @seqs_arr,    $seq );
                push( @headers_arr, $header );

        #				print "\nStoring in index ". ($seq_num - 1) ." $header...\n$seq\n";

                $header_to_index_hash{$header} = ( $seq_num - 1 );
            }
        }

        #		foreach(keys %header_to_index_hash) {
        #			print "header_to_index_hash for $_ is $header_to_index_hash{$_}\n";
        #		}

        print "\n";

        if ( @headers_arr == @seqs_arr ) {

            #			for(my $i=0; $i < @headers_arr; $i++) {
            #				$header_to_index_hash{$headers_arr[$i]} = $i;
            #			}
        }
        else {
            die
"\n\nPROBLEM WITH number of headers not equal to number of sequences\n\n";
        }

        #		print "\nHEADERS:\n";
        #		foreach (@headers_arr) {
        #			print "$_\n";
        #		}
        #
        #		print "\nSEQS:\n";
        #		foreach (@seqs_arr) {
        #			print "$_\n";
        #		}

        ##################################################################################
        # MAIN WORK: CALCULATE NUCLEOTIDE DIVERGENCE dN/dS
        #print "product\tN_sites\tS_sites\tN_diffs\tS_diffs\n";
        my $total_N_sites = 0;
        my $total_S_sites = 0;
        my $total_N_diffs = 0;
        my $total_S_diffs = 0;

        my $between_total_N_sites = 0;
        my $between_total_S_sites = 0;
        my $between_total_N_diffs = 0;
        my $between_total_S_diffs = 0;

        my %min_count2num_alleles_TOTALS;

        #print "\n\n";

        chdir("$OID_USER_DIR\/data\/$family_id");

        # PREPARE FILE (header) of retained positions
        open( OUT_RETAINED_POSNS,
            ">>$OID_USER_DIR/data/$family_id\/retained_aa_positions.txt" );
        print OUT_RETAINED_POSNS "partition\tfamily_positions_retained\n";
        close OUT_RETAINED_POSNS;

        open( WITHIN_GROUP_OUTFILE, ">>within_group_results.txt" );
        print WITHIN_GROUP_OUTFILE "family\tpartition\tN_sites\t"
          . "S_sites\tN_diffs\tS_diffs\t"
          . "dN\tdS\tdN\-dS\tdN\/dS";

        #		if($num_bootstraps > 1) {
        #			print WITHIN_GROUP_OUTFILE "\tSE(dN-dS)\tZ_value\tsignificance\n";
        #		} else {
        print WITHIN_GROUP_OUTFILE "\n";

        #		}

        close WITHIN_GROUP_OUTFILE;

        open( CODON_FILE, ">>within\_group\_codon\_results\.txt" );
        print CODON_FILE
"family\tpartition\tcodon\tvariability\tnum_defined_codons\tcomparisons\t"
          . "N_sites\tS_sites\tN_diffs\tS_diffs\n";
        close CODON_FILE;

        if ($between_group_flag) {

#			open(BETWEEN_CODON_FILE,">>family$family_id\_between\_group\_codon\_results\.txt");
            open( BETWEEN_CODON_FILE, ">>between\_group\_codon\_results\.txt" );
            print BETWEEN_CODON_FILE
"analysis\tfamily\tpartition\tgroup_1\tgroup_2\tcodon\tvariability\t"
              . "num_defined_codons_g1\tnum_defined_codons_g2\tcomparisons\tN_sites\tS_sites\tN_diffs\tS_diffs\n";
            close BETWEEN_CODON_FILE;

   #			open(BETWEEN_OUTFILE,">>family$family_id\_between\_group\_results\.txt");
            open( BETWEEN_OUTFILE, ">>between\_group\_results\.txt" );
            print BETWEEN_OUTFILE
"analysis\tfamily\tpartition\tgroup_1\tgroup_2\tN_sites\tS_sites\tN_diffs\tS_diffs\t"
              . "dN\tdS\t"
              . "dN\-dS\t"
              . "dN\/dS\n";
            close BETWEEN_OUTFILE;

        }

        # FOREACH PARTITION
        my @partitions = keys %{ $family_partition_members_hh{$family_id} };

        foreach my $partition_id (@partitions) {

            print "### Working on partition $partition_id ";
            print
"#########################################################################################\n";

            my $partition_start_time = time;

            if ( $ignore_done_flag == 1 ) {

#				unlink "$OID_USER_DIR/data/$family_id/partition$partition_id\_codon\_results\.txt";
                unlink
"$OID_USER_DIR/data/$family_id/partition$partition_id\_aligned\.faa";
                unlink
"$OID_USER_DIR/data/$family_id/partition$partition_id\_aligned\.cds";

#				unlink "$OID_USER_DIR/data/$family_id/partition$partition_id\_sw_" . $sliding_window_size . "codons\_results.txt";
                unlink
"$OID_USER_DIR/data/$family_id/$partition_id\_bootstrap\_results\.txt";

             # remove lingering <PARTITION#>_polymorphic_codons_arr AND
             # <PARTITION#>_between_polymorphic_codons_arr directory, if present
                if (
                    -d "$OID_USER_DIR\/data\/$family_id\/$partition_id\/$partition_id\_polymorphic_codons_arr"
                  )
                {
                    my @files_to_delete = glob
"$OID_USER_DIR\/data\/$family_id\/$partition_id\/$partition_id\_polymorphic_codons_arr\/*";

                    #					chdir("$OID_USER_DIR\/data\/$family_id");

                    foreach (@files_to_delete) {
                        unlink
"$OID_USER_DIR\/data\/$family_id\/$partition_id\/$partition_id\_polymorphic_codons_arr\/$_";
                    }

                    rmdir
"$OID_USER_DIR\/data\/$family_id\/$partition_id\/$partition_id\_polymorphic_codons_arr";

                }

                if (
                    -d "$OID_USER_DIR\/data\/$family_id\/$partition_id\/$partition_id\_between_polymorphic_codons_arr"
                  )
                {
                    my @files_to_delete = glob
"$OID_USER_DIR\/data\/$family_id\/$partition_id\/$partition_id\_between_polymorphic_codons_arr\/*";

                    #					chdir("$OID_USER_DIR\/data\/$family_id");

                    foreach (@files_to_delete) {
                        unlink
"$OID_USER_DIR\/data\/$family_id\/$partition_id\/$partition_id\_between_polymorphic_codons_arr\/$_";
                    }

                    rmdir
"$OID_USER_DIR\/data\/$family_id\/$partition_id\/$partition_id\_between_polymorphic_codons_arr";
                }

            }

            open( FAMILY_LOG, ">>this_family_log.txt" );
            print FAMILY_LOG "$family_id\t$partition_id\t";

#			open(CODON_FILE,">>partition$partition_id\_codon\_results\.txt");
#			print CODON_FILE "family\tpartition\tcodon\tvariability\tnum_defined_codons\tcomparisons\t".
#				"N_sites\tS_sites\tN_diffs\tS_diffs\n";
#			close CODON_FILE;

            ##############################################################################
      # Read in the partition of amino acid sequences from the FAMILY fasta file
            my $aa_seq = '';
            my %aa_seqs_hash;
            my $aa_header  = '';
            my $aa_seq_num = 0;
            my $aa_last_seq_length;

            open( IN_AA_FASTA, "$this_family_file_name" )
              or die "Could not open file $this_family_file_name\n";

            #EXAMPLE:
            #>Solyc#Solyc01g006920.2.1
            #------------------------------------------------------------
            #------------------------------------------------------------

            print "Recording amino acid data...\n";

            my $include_curr_gene = 0;

            while (<IN_AA_FASTA>) {
                chomp;

                #				print "\nMy inline: $_\n";

                #if(/>([\w\#\.]+])/) {
                if (/>/) {
                    my $match     = $_;
                    my $gene_name = $';    # get rid of '>'

                    #print "match=$match gene_name=$gene_name\n";
                    # match=>Apleuropneumoniae10D13039#appser10_1250
                    # gene_name=Apleuropneumoniae10D13039#appser10_1250

                    # RECORD ONLY IF it's present in the partition
                    if ( $family_partition_members_hh{$family_id}
                        ->{$partition_id}->{$gene_name} == 1 )
                    {

                        $include_curr_gene = 1;

                        if ( $aa_seq_num == 0 ) {
                            $aa_header = $_;
                            $aa_seq_num++;
                        }
                        else {
                            $aa_seqs_hash{$aa_header} = $aa_seq;
                            $aa_header = $_;
                            $aa_seq_num++;

                            my $aa_this_seq_length = length($aa_seq);

                  #print "\nseq $aa_seq_num is of length $aa_this_seq_length\n";
                  #print "\nseq: $aa_seq\n";

                            if (
                                $aa_last_seq_length
                                && (
                                    $aa_last_seq_length != $aa_this_seq_length )
                              )
                            {
                                die
"\n\nDIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
                            }
                            else {
                                $aa_last_seq_length = $aa_this_seq_length;

                                #print "\naa_seq: $aa_seq\n";
                                $aa_seq = '';
                            }
                        }

                    }
                    else
                    {  # ignore this one, because it's not an included partition
                        $include_curr_gene = 0;
                    }

                }
                else {
                    if ( $include_curr_gene == 1 ) {
                        $aa_seq .= $_;
                    }
                }
            }

            close IN_AA_FASTA;

            $aa_seqs_hash{$aa_header} = $aa_seq;

            print "\n";

            ##############################################################################
            # REMOVE amino acid positions with ALL GAPS for this partition;
            # PRINT positions to retained_aa_positions.txt;
            # PRINT alignment to X_aligned.faa
            my %aaPosn_gapCount;
            my $partition_num_seqs = scalar( keys %aa_seqs_hash );

            # Loop through sequences and record number of gaps for each position
            foreach my $curr_header ( keys %aa_seqs_hash ) {
                my $curr_aa_seq = $aa_seqs_hash{$curr_header};

                for (
                    my $index = 0 ;
                    $index < length($curr_aa_seq) ;
                    $index++
                  )
                {
                    my $curr_char = substr( $curr_aa_seq, $index, 1 );

                    if ( $curr_char eq '-' ) {
                        $aaPosn_gapCount{ $index + 1 }++;
                    }
                }
            }

            my @posns_to_remove;
            my $retained_posns_string;
            foreach my $curr_posn ( sort { $a <=> $b } keys %aaPosn_gapCount ) {
                if ( $aaPosn_gapCount{$curr_posn} == $partition_num_seqs )
                {    # all seqs are gap
                     # Store the position for removal
                    push( @posns_to_remove, $curr_posn );
                }
                else {
                    $retained_posns_string .= "$curr_posn,";
                }
            }

            chop($retained_posns_string);

            # Print positions retained to an output file in the family directory
            open( OUT_RETAINED_POSNS,
                ">>$OID_USER_DIR/data/$family_id\/retained_aa_positions.txt" );
            print OUT_RETAINED_POSNS "$partition_id\t$retained_posns_string\n";
            close OUT_RETAINED_POSNS;

        # Loop through sequences, remove the all-gap positions, replace sequence
            my $keep_track_length = 0;
            foreach my $curr_header ( keys %aa_seqs_hash ) {
                my $curr_aa_seq = $aa_seqs_hash{$curr_header};

                my $num_posns_removed = 0;
                for ( my $k = 0 ; $k < scalar(@posns_to_remove) ; $k++ ) {
                    my $posn_to_remove =
                      $posns_to_remove[$k] -
                      $num_posns_removed
                      ; # to account for offset due to already-removed positions
                    my $index_to_remove = $posn_to_remove - 1;
                    substr( $curr_aa_seq, $index_to_remove, 1, '' );
                    $num_posns_removed++;
                }

                if ( $keep_track_length > 0 ) {
                    if ( $keep_track_length != length($curr_aa_seq) ) {
                        die
"\n### DIE: problem removing gaps from amino acid aligments.\n\n";
                    }
                }
                else {    # just first time
                    $keep_track_length = length($curr_aa_seq);
                }

                $aa_seqs_hash{$curr_header} = $curr_aa_seq;

                &print_fasta_file( "partition$partition_id\_aligned.faa",
                    $curr_header, $curr_aa_seq );

            }

            ##############################################################################
         # IMPOSE the amino acid ALIGNMENTS on the relevant nucleotide sequences
            print "Imposing amino acid alignments on nucleotides...\n";
            my @aligned_seqs_arr; # this will be of same length as %aa_seqs_hash
            my @sorted_aa_seqs_hash_keys =
              sort { $a <=> $b } keys %aa_seqs_hash;
            for ( my $i = 0 ; $i < scalar(@sorted_aa_seqs_hash_keys) ; $i++ ) {
                my $this_aln_aa_seq_header = $sorted_aa_seqs_hash_keys[$i];
                my $this_aln_aa_seq = $aa_seqs_hash{$this_aln_aa_seq_header};
                my $this_nt_seq_index;

                #				print "\nThis header is: $this_aln_aa_seq_header\n";

                if ( exists $header_to_index_hash{$this_aln_aa_seq_header} )
                {    # THIS WILL BE ZERO for the ZERO INDEX. Must use EXISTS.
                    $this_nt_seq_index =
                      $header_to_index_hash{$this_aln_aa_seq_header};

                #print "\nWe found $this_aln_aa_seq_header in the .cds files\n";
                }
                else {
                    warn
"\n### WARNING: there is no $this_aln_aa_seq_header in the .cds files\n\n";
                }

                #				print "Imposing alignment on $this_aln_aa_seq_header...\n";

                my $this_nt_seq = $seqs_arr[$this_nt_seq_index];

              #			if($this_aln_aa_seq_header eq '>PmultocidaPm70#NP_246567.1') {
              #				#print "\n\nMy unaligned sequence:\n$this_nt_seq\n";
              #				#print "\n\nMy amino acid sequence:\n$this_aln_aa_seq\n";
              #			}

#			my $this_nt_seq_aligned = &align_codon2aa($this_aln_aa_seq,$this_nt_seq,$this_aln_aa_seq_header);
                my $this_nt_seq_aligned =
                  &align_codon2aa( $this_aln_aa_seq, $this_nt_seq );

                #				print "\n\nMy aligned sequence:\n$this_nt_seq_aligned\n";

                if ($this_nt_seq_aligned)
                {    # if unequal lengths, etc., nothing returned
                     # Print (append) the aligned nucleotide sequence to a file in the family's folder
##					&print_fasta_file("aligned.cds",$this_aln_aa_seq_header,$this_nt_seq_aligned);
                    &print_fasta_file( "partition$partition_id\_aligned.cds",
                        $this_aln_aa_seq_header, $this_nt_seq_aligned );

                    push( @aligned_seqs_arr, $this_nt_seq_aligned );
                }
                else {
                    warn
"\n### WARNING: The protein sequence $this_aln_aa_seq_header does not contain "
                      . "the expected number of amino acids based on the length of nucleotide "
                      . "sequence. Sequence excluded from dN/dS analysis.\n\n";
                }
            }

            ##############################################################################
            # FOR EACH ALIGNED SEQUENCE: perform dN/dS
            #			print "\nAnalyzing $family_id:\n";

            my @codonum_codon_aa;
            my $last_num_codons;

            # ADD: minimum freq and count for pairwise deletion
            my @codonIndex2count;

            # STORE all CODONS here
          STORE_CODONS:
            for (
                my $seq_id = 0 ;
                $seq_id < scalar(@sorted_aa_seqs_hash_keys) ;
                $seq_id++
              )
            {    # for each sequence (in partition?)

#print "\nValue of sorted_aa_seqs_hash_keys [ seq_id==$seq_id ]: $sorted_aa_seqs_hash_keys[$seq_id]\n";
# For example:
# Value of sorted_aa_seqs_hash_keys [ seq_id==19 ]: >Msuccinici#YP_088511.1

                my $product_seq    = $aligned_seqs_arr[$seq_id];
                my $product_length = length($product_seq);

                # Make sure it's a complete set of codons
                if ( ( $product_length % 3 ) != 0 ) {

#die "\n\nDIE: A sequence in $family_id is not a multiple of 3 (complete codon set). TERMINATED.\n\n";
                    warn
"\n### WARNING: A sequence in family $family_id partition $partition_id is not a multiple of 3 (complete codon set). "
                      . "Sequence ignored in dN/dS analysis.\n\n";
                    next STORE_CODONS;
                }

                # Build all codons
                my $num_codons = $product_length / 3;

                if (   ( !$last_num_codons )
                    && ( $last_num_codons ne '' )
                    && ( $num_codons != $last_num_codons ) )
                {
#die "\n\nDIE: In $family_id, there are sequences of different length ($last_num_codons and $num_codons). TERMINATED.\n\n";
#				print "\n\nWARNING: In $family_id, there are sequences of different length ($last_num_codons and $num_codons).".
#					"Sequence ignored in dN/dS analysis.\n\n";
                    next STORE_CODONS;
                }
                else {
                    $last_num_codons = $num_codons;
                }

                my $codon_index = 0;    # populate @codonum_codon_aa
                for ( my $i = 1 ; $i <= $num_codons ; $i++ ) {  # for each codon
                    my $array_index = $i - 1;

                    #my $codon = substr($sequence,$codon_index,3);
                    # More efficient to call substr at beginning and gobble
                    my $codon = substr( $product_seq, 0, 3, "" );

                    #print "$codon ";

                    # Note that some may contain N's, gaps (-)
                    push( @{ $codonum_codon_aa[$array_index] }, $codon );

                    #$codon_index+=3;

                    # ADDED if the codon doesn't contain a gap or N
                    if ( !( $codon =~ /-/ || $codon =~ /N/ || $codon eq '' ) ) {
                        $codonIndex2count[$array_index]++;
                    }
                }

                #print "\n";
            }  # STORE_CODONS; end sequences loop; end STORING codon information

            ##############################################################################
            # ANALYZE codons here; all pairwise comparisons

            if ( !@codonum_codon_aa ) {
                print "\nFamily $family_id partition $partition_id ignored\n";

                #next FAMILY; # can't do with parallel fork

                open( IGNOREFILE,
                    ">>$OID_USER_DIR/data/$family_id/.snpgenie_ignore" );
                print IGNOREFILE "done";
                close IGNOREFILE;

            }
            else {

                my $num_seqs = scalar @{ $codonum_codon_aa[0] };
                my @num_seqs_defined;

                #print "\nNum seqs: $num_seqs\n";
                my $num_codons_in_product = scalar @codonum_codon_aa;
                my @polymorphic_codons_arr
                  ;    # same indices as codons; value 0 if invariant, 1 if poly

                print FAMILY_LOG "$num_codons_in_product\t";

     # QUICKLY CHECK IF IT'S POLYMORPHIC to save time; FOR EACH CODON IN PRODUCT
                for (
                    my $codon_index = 0 ;
                    $codon_index < $num_codons_in_product ;
                    $codon_index++
                  )
                {

               #print "Codon $codon_index\n";
               #my $codon_to_print = $groups_codons_aa[$i]->[$codon_index]->[1];
               #$| = 1;
               #print "$codon_to_print";
                    my $between_s_comparisons       = 0;
                    my $this_codon_poly             = 0;
                    my $this_codon_num_seqs_defined = $num_seqs;

                    # Check num defined sequences
                    for (
                        my $si_seq_index = 0 ;
                        $si_seq_index < $num_seqs ;
                        $si_seq_index++
                      )
                    {    # for each si sequence at this codon
                        my $codon_si =
                          $codonum_codon_aa[$codon_index]->[$si_seq_index];

                        # Undefined?
                        if ( $codon_si =~ "[^ACGT]" ) {
                            $this_codon_num_seqs_defined -=
                              1;    # will this really work here? COMEBACK NOW
                        }
                    }
                    push( @num_seqs_defined, $this_codon_num_seqs_defined );

                    # Check polymorphic
                  OUTER:
                    for (
                        my $si_seq_index = 0 ;
                        $si_seq_index < $num_seqs ;
                        $si_seq_index++
                      )
                    {    # for each si sequence at this codon
                        my $codon_si =
                          $codonum_codon_aa[$codon_index]->[$si_seq_index];

                        # Codon to compare against
                        if ( !( $codon_si =~ 'N' ) && !( $codon_si =~ '-' ) ) {

                          INNER:
                            for (
                                my $sj_seq_index = $si_seq_index + 1 ;
                                $sj_seq_index < $num_seqs ;
                                $sj_seq_index++
                              )
                            {    # for each sj sequence at this codon
                                 #if ($codon_index == 0) { print "pw comparison\n"; }
                                my $codon_sj = $codonum_codon_aa[$codon_index]
                                  ->[$sj_seq_index];

                            #print "Comparing codons $codon_si and $codon_sj\n";
                                if (   $codon_si ne $codon_sj
                                    && !( $codon_sj =~ 'N' )
                                    && !( $codon_sj =~ '-' ) )
                                {
                                    $this_codon_poly = 1;

#							if($family_id == 80 && $codon_index >= 546) {
#								print "\nFamily $family_id codon_index $codon_index: polymorphic because $codon_si versus $codon_sj\n";
#							}

                                    last OUTER
                                      ; # gotta break of out TWO loops here and go to next codon in product
                                }
                            }    # INNER
                        }    #else {
                         #	$this_codon_num_seqs_defined -= 1; # will this really work here? COMEBACK NOW
                         #}
                    }    # OUTER

                    # Save results
                    push( @polymorphic_codons_arr, $this_codon_poly );

                  #$num_seqs_defined[$poly_file] = $this_codon_num_seqs_defined;

                }    # done checking if polymorphic

                print "\nAnalyzing polymorphic codons...\n";

# STORE DATA HERE (this partition) for later between-group and sliding window analyses
                my %data_codonnum_spp_codon_hh;

                my $product_N_sites_sum = 0;
                my $product_S_sites_sum = 0;
                my $product_N_diffs_sum = 0;
                my $product_S_diffs_sum = 0;

          #				open(CODON_FILE,">>partition$partition_id\_codon\_results\.txt");
                open( CODON_FILE, ">>within\_group\_codon\_results\.txt" );

# ADD: COMEBACK: make a number of sites subroutine for the whole codon, not just the sites
# OR, instead of calling the subroutine every time, just COUNT the comparisons!

                # POLYMORPHIC SITE ANALYSIS
                my %poly_site_data_hh;
                my %codon_to_results_hh;    # for bootstraps
                my $output_buffer = '';

              COMPL_DEL:
                for (
                    my $codon_index = 0 ;
                    $codon_index < $num_codons_in_product ;
                    $codon_index++
                  )
                {                           # for each codon in product
                    my $codon_num = $codon_index + 1;

                    $output_buffer .=
                      "$family_id\t$partition_id\t" . $codon_num . " \t";

           #print CODON_FILE "$family_id\t$partition_id\t" . $codon_num . " \t";

                    # Between-group storage
###					my %species_multi_hits;

                    if ($between_group_flag) {
###						my %seen_species_hash;

                        for (
                            my $seq_index = 0 ;
                            $seq_index < $num_seqs ;
                            $seq_index++
                          )
                        {    # for each si sequence at this codon
                            my $codon_si =
                              $codonum_codon_aa[$codon_index]->[$seq_index];

                            my $species_i;

                           # HERE IT IS. We're assuming ONLY ONE SPECIES MEMBER.
                            if ( $sorted_aa_seqs_hash_keys[$seq_index] =~
                                />(\w+)#/ )
                            {
                                $species_i = $1;

###								if($seen_species_hash{$species_i} == 1) {
###									$seen_species_hash{$species_i} = 1;
###								} else {
###									$species_multi_hits{$species_i} = 1;
###								}

                                $data_codonnum_spp_codon_hh{$codon_num}
                                  ->{$species_i} = $codon_si
                                  ;    # STORING ALL — even --- or NNN or ''
                            }
                            else {
                                die "\n### DIE: species name not found.\n\n";
                            }
                        }

                    }    # end between-group storage

###					foreach(sort keys %species_multi_hits) {
###						warn "### WARNING: there are multiple members of species $_ " .
###							"in family $family_id partition $partition_id.\n" .
###							"### Please fix your data. There should be one representative ortholog for each species in a partition.\n" .
###							"### One gene will be arbitrary chosen to proceed.\n\n";
###					}

                    # ADD: minimum frequency/codon in pairwise deletion method
                    if ( $polymorphic_codons_arr[$codon_index] )
                    {    # if codon is polymorphic, value is 1

                        my $this_actual_count = $codonIndex2count[$codon_index];

                        #				$min_freq_in_alignment;
                        #				$min_count_in_alignment;

                        # Make sure we make the minimum count/frequency
                        #						if($this_actual_count >= $this_min_count) {
                        #print "\nGot here 5\n";
                        $output_buffer .= "polymorphic\t$this_actual_count\t";

                        #print CODON_FILE "polymorphic\t$this_actual_count\t";

                        my %comps_hh;

                        for (
                            my $si_seq_index = 0 ;
                            $si_seq_index < $num_seqs ;
                            $si_seq_index++
                          )
                        {    # for each si sequence at this codon
                            my $codon_si =
                              $codonum_codon_aa[$codon_index]->[$si_seq_index];

                            # Sequence (codon) to compare against
                            if (   !( $codon_si =~ 'N' )
                                && !( $codon_si =~ '-' ) )
                            {

                                for (
                                    my $sj_seq_index = $si_seq_index + 1 ;
                                    $sj_seq_index < $num_seqs ;
                                    $sj_seq_index++
                                  )
                                {    # for each sj sequence at this codon
                                    my $codon_sj =
                                      $codonum_codon_aa[$codon_index]
                                      ->[$sj_seq_index];

                                    if (   !( $codon_sj =~ 'N' )
                                        && !( $codon_sj =~ '-' ) )
                                    { # they don't have to be the same; syn still stored here
                                        $comps_hh{$codon_si}->{$codon_sj} += 1;

                                    }
                                }
                            }
                        }

                        # SUM UP STUFF HERE
                        my $between_s_comparisons = 0;
                        my $sum_N_sites           = 0;
                        my $sum_S_sites           = 0;
                        my $sum_N_diffs           = 0;
                        my $sum_S_diffs           = 0;

                        foreach my $codon_si ( keys %comps_hh ) {
                            foreach
                              my $codon_sj ( keys %{ $comps_hh{$codon_si} } )
                            {

                                # ADD OR REVISE
                                if (   $codon_si =~ /-/
                                    || $codon_sj =~ /-/
                                    || $codon_si =~ /N/
                                    || $codon_sj =~ /N/ )
                                {    # these will skip all additions
                                    print "\nWe got here\n";
                                    $output_buffer .=
                                      "COMPLETE_DELETION\tNA\tNA\tNA\tNA\n";

                        #print CODON_FILE "COMPLETE_DELETION\tNA\tNA\tNA\tNA\n";
                                    next COMPL_DEL;
                                }
                                else {

                                    my $weight =
                                      $comps_hh{$codon_si}->{$codon_sj};

                                    $poly_site_data_hh{$codon_index}
                                      ->{$codon_si} += $weight;
                                    $poly_site_data_hh{$codon_index}
                                      ->{$codon_sj} += $weight;

                                    $output_buffer .= "($codon_si\-$codon_sj)";

                                    #print CODON_FILE "($codon_si\-$codon_sj)";

                                    # Sites codon gi
                                    my @codon_si_sites_1_arr =
                                      &get_number_of_sites( $codon_si, 1 );
                                    my @codon_si_sites_2_arr =
                                      &get_number_of_sites( $codon_si, 2 );
                                    my @codon_si_sites_3_arr =
                                      &get_number_of_sites( $codon_si, 3 );

                                    my $codon_si_N_sites =
                                      ( $codon_si_sites_1_arr[0] +
                                          $codon_si_sites_2_arr[0] +
                                          $codon_si_sites_3_arr[0] );
                                    my $codon_si_S_sites =
                                      ( $codon_si_sites_1_arr[1] +
                                          $codon_si_sites_2_arr[1] +
                                          $codon_si_sites_3_arr[1] );

                                    # Site codon gj
                                    my @codon_sj_sites_1_arr =
                                      &get_number_of_sites( $codon_sj, 1 );
                                    my @codon_sj_sites_2_arr =
                                      &get_number_of_sites( $codon_sj, 2 );
                                    my @codon_sj_sites_3_arr =
                                      &get_number_of_sites( $codon_sj, 3 );

                                    my $codon_sj_N_sites =
                                      ( $codon_sj_sites_1_arr[0] +
                                          $codon_sj_sites_2_arr[0] +
                                          $codon_sj_sites_3_arr[0] );
                                    my $codon_sj_S_sites =
                                      ( $codon_sj_sites_1_arr[1] +
                                          $codon_sj_sites_2_arr[1] +
                                          $codon_sj_sites_3_arr[1] );

                         #print "\nComparing codons $codon_si and $codon_sj:\n";

                                    my $mean_comp_N_sites =
                                      ( $codon_si_N_sites + $codon_sj_N_sites )
                                      / 2;
                                    my $mean_comp_S_sites =
                                      ( $codon_si_S_sites + $codon_sj_S_sites )
                                      / 2;

                                    $sum_N_sites +=
                                      $weight * $mean_comp_N_sites;
                                    $sum_S_sites +=
                                      $weight * $mean_comp_S_sites;

                                    # Differences
                                    my $N_diffs = 0;
                                    my $S_diffs = 0;

                                    if ( $codon_si ne $codon_sj ) {
                                        my @diffs_arr =
                                          &return_avg_diffs( $codon_si,
                                            $codon_sj );
                                        $N_diffs = $diffs_arr[0];
                                        $S_diffs = $diffs_arr[1];
                                    }

                                    $sum_N_diffs += $weight * $N_diffs;
                                    $sum_S_diffs += $weight * $S_diffs;

                                    $between_s_comparisons += $weight;
                                }
                            }
                        }

                        $output_buffer .= "\t";

                        #print CODON_FILE "\t";

                 #print "\nfor a total of $between_s_comparisons comparisons\n";

                        my $mean_N_sites =
                          ( $sum_N_sites / $between_s_comparisons );
                        my $mean_S_sites =
                          ( $sum_S_sites / $between_s_comparisons );
                        my $mean_N_diffs =
                          ( $sum_N_diffs / $between_s_comparisons );
                        my $mean_S_diffs =
                          ( $sum_S_diffs / $between_s_comparisons );

                        $output_buffer .=
"$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs\n";

#print CODON_FILE "$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs\n";

                        $product_N_sites_sum += $mean_N_sites;
                        $product_S_sites_sum += $mean_S_sites;
                        $product_N_diffs_sum += $mean_N_diffs;
                        $product_S_diffs_sum += $mean_S_diffs;

            #print ".. and it's all done.\n";
            #						} # else we didn't make the minimum count/frequency threshold

                    }
                    else
                    { # not polymorphic; just use first codon from the first group because it's conserved
                        my $this_actual_count = $codonIndex2count[$codon_index];

                        my $conserved_codon = '';

                      FIND_CONSERVED_CODON:
                        for ( my $k = 0 ; $k < $num_seqs ; $k++ )
                        {    # for each codon in product
                            my $curr_codon_rep =
                              $codonum_codon_aa[$codon_index]->[$k]
                              ;    # grab first codon THAT DOESN'T CONTAIN N
                            if (   $curr_codon_rep =~ "N"
                                || $curr_codon_rep =~ "-"
                                || $curr_codon_rep eq '' )
                            {
                                next FIND_CONSERVED_CODON;
                            }
                            else {
                                $conserved_codon = $curr_codon_rep;
                                last FIND_CONSERVED_CODON;
                            }
                        }

                        my @codon_sites_1_arr =
                          &get_number_of_sites( $conserved_codon, 1 );
                        my @codon_sites_2_arr =
                          &get_number_of_sites( $conserved_codon, 2 );
                        my @codon_sites_3_arr =
                          &get_number_of_sites( $conserved_codon, 3 );

                        my $codon_N_sites =
                          ( $codon_sites_1_arr[0] +
                              $codon_sites_2_arr[0] +
                              $codon_sites_3_arr[0] );
                        my $codon_S_sites =
                          ( $codon_sites_1_arr[1] +
                              $codon_sites_2_arr[1] +
                              $codon_sites_3_arr[1] );

                        if (   $this_actual_count == 0
                            || $this_actual_count eq ''
                            || !defined($this_actual_count) )
                        {
                            $this_actual_count = 0;
                        }

                        $output_buffer .=
                            "conserved\t$this_actual_count\t$conserved_codon\t"
                          . "$codon_N_sites\t$codon_S_sites\t0\t0\n";

          #print CODON_FILE "conserved\t$this_actual_count\t$conserved_codon\t".
          #	"$codon_N_sites\t$codon_S_sites\t0\t0\n";

                        $product_N_sites_sum += $codon_N_sites;
                        $product_S_sites_sum += $codon_S_sites;

                    }    # end poly/not poly check

                    if ( length($output_buffer) > 50000 ) {
                        print CODON_FILE "$output_buffer";
                        $output_buffer = '';
                    }

                }    # end all codons in product

                # CLEAR BUFFER
                print CODON_FILE "$output_buffer";
                $output_buffer = '';
                close CODON_FILE;

#				# Check contents for between-group analysis
#				for(my $codon_index = 0; $codon_index < $num_codons_in_product; $codon_index++) { # for each codon in product
#
#					my $codon_num = $codon_index+1;
#
#					print "\nCODON: $codon_num\n";
#
#					foreach(keys %{$data_codonnum_spp_codon_hh{$codon_num}}) {
#
#						my $species = $_;
#						my $codon = $data_codonnum_spp_codon_hh{$codon_num}->{$species};
#						print "species=$species codon=$codon\n";
#
#					}
#				}

                ##########################################################################
                # OPTIONAL between-group analysis
                if ($between_group_flag) {

                    if ( -f "$OID_USER_DIR/between_group.config" ) {

                        print
"\n################### Performing between-group analyses...\n";

                        open( BETWEEN_GROUP_CONFIG,
                            "$OID_USER_DIR/between_group.config" )
                          or die
"Could not open file $OID_USER_DIR/between_group.config\n";

                        while (<BETWEEN_GROUP_CONFIG>)
                        {    # for each analysis requested
                            chomp;

                            my $this_line = $_;

                            my @this_line_arr = split( /\s+/, $this_line );

                            my $analysis_name = $this_line_arr[0];

      # Get species for each group.
      # Remember, input files contain EVERY SPECIES POSSIBLE...
      # we want to cull them down to only the SPECIES PRESENT IN THIS PARTITION.
                            my @group_1_spp;
                            my @group_2_spp;

                            if ( $this_line_arr[1] =~ /\(([\w\,]+)\)/ ) {
                                my $species_string = $1;
                                my @group_1_possible_spp =
                                  split( /,/, $species_string );

                                my %added_species;

                                foreach (@group_1_possible_spp) {

                                    if (   $family_species_h{$_} == 1
                                        && $added_species{$_} != 1 )
                                    {
                                        $added_species{$_} = 1;
                                        push( @group_1_spp, $_ );
                                    }
                                }
                            }

                            if ( $this_line_arr[2] =~ /\(([\w\,]+)\)/ ) {
                                my $species_string = $1;
                                my @group_2_possible_spp =
                                  split( /,/, $species_string );

                                my %added_species;

                                foreach (@group_2_possible_spp) {

                                    if (   $family_species_h{$_} == 1
                                        && $added_species{$_} != 1 )
                                    {
                                        $added_species{$_} = 1;
                                        push( @group_2_spp, $_ );
                                    }
                                }
                            }

                    # We are comparing two groups: @group_1_spp and @group_2_spp
                            print
"\nAnalysis: $analysis_name for $num_codons_in_product codons in partition $partition_id\...\n";

                            my $group_name_i = "(@group_1_spp\)";
                            my $group_name_j = "(@group_2_spp\)";

#print "\nmy group name i is: " . $fasta_files[0] . "\n";
#print "\nmy group name i is: $group_name_i\nmy group name j is: $group_name_j\n";

                        # STORE THE SEQUENCES FOR THIS PRODUCT FROM GROUP J ONLY
                        # the fasta file corresponding to group j is:

                            # Updated here:
                            my $num_seqs_per_codon_groupi = scalar @group_1_spp;
                            my $num_seqs_per_codon_groupj = scalar @group_2_spp;

# Quickly see if it's polymorphic between-group to save time. It's faster just to check, EVEN IF not polymorphic
                            my @between_polymorphic_codons_arr;
                            my %codonIndex_conserved;

                            print
                              "\nDetermining which codons are polymorphic...\n";

                # DETERMINE WHICH CODONS ARE POLYMORPHIC, **BETWEEN** the GROUPS
                            for (
                                my $codon_index = 0 ;
                                $codon_index < $num_codons_in_product ;
                                $codon_index++
                              )
                            {    # for each codon in product

                                my $codon_num       = $codon_index + 1;
                                my $this_codon_poly = 0;

                              OUTER:
                                for (
                                    my $gi_spp_index = 0 ;
                                    $gi_spp_index < $num_seqs_per_codon_groupi ;
                                    $gi_spp_index++
                                  )
                                {    # for each gi sequence at this codon

                                    # Codon to compare against
                                  INNER:
                                    for (
                                        my $gj_spp_index = 0 ;
                                        $gj_spp_index <
                                        $num_seqs_per_codon_groupj ;
                                        $gj_spp_index++
                                      )
                                    {    # for each gj sequence at this codon

                                        my $group_1_species =
                                          $group_1_spp[$gi_spp_index];
                                        my $group_2_species =
                                          $group_2_spp[$gj_spp_index];

                                        my $codon_gi =
                                          $data_codonnum_spp_codon_hh{
                                            $codon_num}->{$group_1_species};
                                        my $codon_gj =
                                          $data_codonnum_spp_codon_hh{
                                            $codon_num}->{$group_2_species};

                            #print "Comparing codons $codon_gi and $codon_gj\n";
                                        if (
                                            ( $codon_gi ne $codon_gj )
                                            && !(
                                                   $codon_gi =~ 'N'
                                                || $codon_gi =~ '-'
                                                || $codon_gi eq ''
                                                || $codon_gj =~ 'N'
                                                || $codon_gj =~ '-'
                                                || $codon_gj eq ''
                                            )
                                          )
                                        {
                                            $this_codon_poly = 1;

#print "\nbetween-group codon index $codon_index is polymorphic\n".
#	"species $group_1_species vs. $group_2_species\ncodon $codon_gi vs. $codon_gj\n";
                                            last OUTER
                                              ; # gotta break of out TWO loops here and go to next codon in product
                                        }    #else {
                                         #	print "\nThis is conserved: codon $codon_num ($codon_gi)/($codon_gj)\n";
                                         #}
                                    }    # INNER

                                }    # OUTER

                                # Save results
                                push( @between_polymorphic_codons_arr,
                                    $this_codon_poly );

                            }

                            print
"\nDone recording polymorphism in partition $partition_id between the groups:\nGROUP 1: $group_name_i\nGROUP 2: $group_name_j\n";

                            #print "It is: @between_polymorphic_codons_arr\n\n";

                            print
"\nAnalyzing polymorphic codons in partition $partition_id\...\n";

          # calculate number of total codons so that we can track % progress?
          # INTRO: skip if not poly, just do values for the first codon observed

                            print "\nAnalyzing each codon...\n";

                            my %codonum_group_numDefinedCodons_hh;

                            # STORE FOR SLIDING WINDOWS AND BOOTSTRAPPING
                            my %codon_data_hh
                              ; # Will store codon->N_diffs_codon/S_diffs_codon/...
                            my @codon_lines_arr;

                          BETWEEN_PW_DEL:
                            for (
                                my $codon_index = 0 ;
                                $codon_index < $num_codons_in_product ;
                                $codon_index++
                              )
                            {    # for each codon in product

                                my $codon_num         = $codon_index + 1;
                                my $this_codon_output = '';

                               #print "codon $codon_num at time " . time . "\n";

#print CODON_FILE "$partition_id\t$i\t$j\tcodon_" . $codon_num . " \t";
#				print CODON_FILE "$partition_id\t$group_name_i\t$group_name_j\tcodon_" . $codon_num . " \t";

#"analysis\tfamily\tpartition\tgroup_1\tgroup_2\tcodon\tvariability\tcomparisons\tN_sites\tS_sites\tN_diffs\tS_diffs\n";

                                $this_codon_output .=
"$analysis_name\t$family_id\t$partition_id\t$group_name_i\t$group_name_j\t$codon_num\t";

# We have stored:
# $data_codonnum_spp_codon_hh{$codon_num}->{$species_i} = $codon_si; # STORING ALL — even --- or NNN or ''
# This is everything after the > and before the # (gene name)

# COMEBACK: use ! '' || [^ACGT] and test to make sure
# COMEBACK: also, both here and earlier for defined sites, let's jjust go through and see if it mismatches the LAST DEFINED one. Don't need all pairs!
# Count number of defined codons in group i
                                for (
                                    my $gi_spp_index = 0 ;
                                    $gi_spp_index < $num_seqs_per_codon_groupi ;
                                    $gi_spp_index++
                                  )
                                {    # for each gj sequence at this codon
                                    my $group_1_species =
                                      $group_1_spp[$gi_spp_index];

                                    my $this_species_codon =
                                      $data_codonnum_spp_codon_hh{$codon_num}
                                      ->{$group_1_species};

                                    if (
                                        !(
                                               $this_species_codon =~ /N/
                                            || $this_species_codon =~ /-/
                                            || $this_species_codon eq ''
                                        )
                                      )
                                    {
                                        $codonum_group_numDefinedCodons_hh{
                                            $codon_num}->{group1}++;
                                    }
                                }

                                # Count number of defined codons in group j
                                for (
                                    my $gj_spp_index = 0 ;
                                    $gj_spp_index < $num_seqs_per_codon_groupj ;
                                    $gj_spp_index++
                                  )
                                {    # for each gj sequence at this codon
                                    my $group_2_species =
                                      $group_2_spp[$gj_spp_index];

                                    my $this_species_codon =
                                      $data_codonnum_spp_codon_hh{$codon_num}
                                      ->{$group_2_species};

                                    if (
                                        !(
                                               $this_species_codon =~ /N/
                                            || $this_species_codon =~ /-/
                                            || $this_species_codon eq ''
                                        )
                                      )
                                    {
                                        $codonum_group_numDefinedCodons_hh{
                                            $codon_num}->{group2}++;
                                    }
                                }

                                #my $curr_comp_num = 0;

                 # ONLY IF POLYMORPHIC between-group (WILL THERE BE DIFFERENCES)
                                if (
                                    $between_polymorphic_codons_arr[
                                    $codon_index]
                                  )
                                {    # polymorphic codon
                                     #					print CODON_FILE "polymorphic\t";
                                    $this_codon_output .=
                                      "polymorphic\t"
                                      . $codonum_group_numDefinedCodons_hh{
                                        $codon_num}->{group1}
                                      . "\t"
                                      . $codonum_group_numDefinedCodons_hh{
                                        $codon_num}->{group2}
                                      . "\t";

                                    my %between_comps_hh;

                                    # Sequence (codon) from group 1
                                    for (
                                        my $gi_spp_index = 0 ;
                                        $gi_spp_index <
                                        $num_seqs_per_codon_groupi ;
                                        $gi_spp_index++
                                      )
                                    {    # for each gi sequence at this codon

                                  # Sequence (codon) to compare against, group 2
                                        for (
                                            my $gj_spp_index = 0 ;
                                            $gj_spp_index <
                                            $num_seqs_per_codon_groupj ;
                                            $gj_spp_index++
                                          )
                                        {   # for each gj sequence at this codon
                                            my $group_1_species =
                                              $group_1_spp[$gi_spp_index];
                                            my $group_2_species =
                                              $group_2_spp[$gj_spp_index];

                                            my $codon_gi =
                                              $data_codonnum_spp_codon_hh{
                                                $codon_num}->{$group_1_species};
                                            my $codon_gj =
                                              $data_codonnum_spp_codon_hh{
                                                $codon_num}->{$group_2_species};

#print "\ngroup_1_species is $group_1_species group_2_species is $group_2_species\ncodon_gi is $codon_gi codon_gj is $codon_gj\n";

                                            if (
                                                !(
                                                       $codon_gj =~ 'N'
                                                    || $codon_gj =~ '-'
                                                    || $codon_gj eq ''
                                                    || $codon_gi =~ 'N'
                                                    || $codon_gi =~ '-'
                                                    || $codon_gi eq ''
                                                )
                                              )
                                            { # they don't have to be the same; syn still stored here
                                                 #if(($codon_gi ne '') && ($codon_gj ne '')) {
                                                 #print "we stored $codon_gi and $codon_gj\n";
                                                $between_comps_hh{$codon_gi}
                                                  ->{$codon_gj} += 1;
                                            }
                                        }
                                    }

                                    # ^ do this_actual_count here above?

                                    # SUM UP STUFF HERE
                                    my $between_g_comparisons = 0;
                                    my $between_sum_N_sites   = 0;
                                    my $between_sum_S_sites   = 0;
                                    my $between_sum_N_diffs   = 0;
                                    my $between_sum_S_diffs   = 0;

                                  #print"\nBetween groups $i and $j we compare";

                                    foreach
                                      my $codon_gi ( keys %between_comps_hh )
                                    {

                                        foreach my $codon_gj (
                                            keys
                                            %{ $between_comps_hh{$codon_gi} } )
                                        {

                                            # What's going on here? COMEBACK
                                            if (
                                                (
                                                       $codon_gi =~ /-/
                                                    || $codon_gi =~ /N/
                                                )
                                                && (   $codon_gj =~ /-/
                                                    || $codon_gj =~ /N/ )
                                              )
                                            {    # these will skip all additions
                                                 #print "\nboth contained - or N\n";
                                                 #print CODON_FILE "COMPLETE_DELETION\tNA\tNA\tNA\tNA\n";
                                                 #$this_codon_output .= "PAIRWISE_DELETION\tNA\tNA\tNA\tNA\n"; # this was pointless because it's a private string
                                                next BETWEEN_PW_DEL;
                                            }
                                            else {
                                                my $weight =
                                                  $between_comps_hh{$codon_gi}
                                                  ->{$codon_gj};

                                                ### NEW POLYMORPHIC SITE ANALYSIS ADDED
              #								$poly_site_data_hh{$codon_index}->{$codon_si} += $weight;
              #								$poly_site_data_hh{$codon_index}->{$codon_sj} += $weight;
                                                ###

                             #								print CODON_FILE "($codon_gi\-$codon_gj)";
                                                $this_codon_output .=
                                                  "($codon_gi\-$codon_gj)";

                                                # Sites codon gi
                                                my @codon_gi_sites_1_arr =
                                                  &get_number_of_sites(
                                                    $codon_gi, 1 );
                                                my @codon_gi_sites_2_arr =
                                                  &get_number_of_sites(
                                                    $codon_gi, 2 );
                                                my @codon_gi_sites_3_arr =
                                                  &get_number_of_sites(
                                                    $codon_gi, 3 );

                                                my $codon_gi_N_sites =
                                                  ( $codon_gi_sites_1_arr[0] +
                                                        $codon_gi_sites_2_arr[0]
                                                      + $codon_gi_sites_3_arr[0]
                                                  );
                                                my $codon_gi_S_sites =
                                                  ( $codon_gi_sites_1_arr[1] +
                                                        $codon_gi_sites_2_arr[1]
                                                      + $codon_gi_sites_3_arr[1]
                                                  );

                                                # Site codon gj
                                                my @codon_gj_sites_1_arr =
                                                  &get_number_of_sites(
                                                    $codon_gj, 1 );
                                                my @codon_gj_sites_2_arr =
                                                  &get_number_of_sites(
                                                    $codon_gj, 2 );
                                                my @codon_gj_sites_3_arr =
                                                  &get_number_of_sites(
                                                    $codon_gj, 3 );

                                                my $codon_gj_N_sites =
                                                  ( $codon_gj_sites_1_arr[0] +
                                                        $codon_gj_sites_2_arr[0]
                                                      + $codon_gj_sites_3_arr[0]
                                                  );
                                                my $codon_gj_S_sites =
                                                  ( $codon_gj_sites_1_arr[1] +
                                                        $codon_gj_sites_2_arr[1]
                                                      + $codon_gj_sites_3_arr[1]
                                                  );

                         #print "\nComparing codons $codon_gi and $codon_gj:\n";

                                                my $mean_comp_N_sites =
                                                  ( $codon_gi_N_sites +
                                                      $codon_gj_N_sites ) / 2;
                                                my $mean_comp_S_sites =
                                                  ( $codon_gi_S_sites +
                                                      $codon_gj_S_sites ) / 2;

                                                $between_sum_N_sites +=
                                                  $weight * $mean_comp_N_sites;
                                                $between_sum_S_sites +=
                                                  $weight * $mean_comp_S_sites;

                                                # Differences
                                                my $N_diffs = 0;
                                                my $S_diffs = 0;

                                                if ( $codon_gi ne $codon_gj ) {
                                                    my @diffs_arr =
                                                      &return_avg_diffs(
                                                        $codon_gi, $codon_gj );
                                                    $N_diffs = $diffs_arr[0];
                                                    $S_diffs = $diffs_arr[1];
                                                }

                                                $between_sum_N_diffs +=
                                                  $weight * $N_diffs;
                                                $between_sum_S_diffs +=
                                                  $weight * $S_diffs;

                                                $between_g_comparisons +=
                                                  $weight;

#								$this_codon_data .= "$weight\t$mean_comp_N_sites\t$mean_comp_S_sites\t" .
#									"$N_diffs\t$S_diffs";
                                            }
                                        }
                                    }

                                    #					print CODON_FILE "\t";
                                    $this_codon_output .= "\t";

                 #print "\nfor a total of $between_g_comparisons comparisons\n";

                                    my $mean_N_sites = 0;
                                    my $mean_S_sites = 0;
                                    my $mean_N_diffs = 0;
                                    my $mean_S_diffs = 0;

                                    if ( $between_g_comparisons > 0 ) {
                                        $mean_N_sites = ( $between_sum_N_sites /
                                              $between_g_comparisons );
                                        $mean_S_sites = ( $between_sum_S_sites /
                                              $between_g_comparisons );
                                        $mean_N_diffs = ( $between_sum_N_diffs /
                                              $between_g_comparisons );
                                        $mean_S_diffs = ( $between_sum_S_diffs /
                                              $between_g_comparisons );
                                    }

#					print CODON_FILE "$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs\n";
                                    $this_codon_output .=
"$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs";

                                    # SLIDING WINDOW & BOOTSTRAP
                                    $codon_data_hh{$codon_num}->{N_sites} =
                                      $mean_N_sites;
                                    $codon_data_hh{$codon_num}->{S_sites} =
                                      $mean_S_sites;
                                    $codon_data_hh{$codon_num}->{N_diffs} =
                                      $mean_N_diffs;
                                    $codon_data_hh{$codon_num}->{S_diffs} =
                                      $mean_S_diffs;

                                    ###
                                    #					$product_N_sites_sum += $mean_N_sites;
                                    #					$product_S_sites_sum += $mean_S_sites;
                                    #					$product_N_diffs_sum += $mean_N_diffs;
                                    #					$product_S_diffs_sum += $mean_S_diffs;
                                    ###

                                    #print ".. and it's all done.\n";

                                }
                                else
                                { # not polymorphic; just use first DEFINED codon from the first group because it's conserved
                                     #my $conserved_codon = $groups_codons_aa[$i]->[$codon_index]->[0]; # grab first codon10
                                     # BUT, if the first group is all gaps...?

                                    my $conserved_codon = '';

#									if($codon_num == 2 && $partition_id == 1) {
#										print "\nAccording to my records, there are $num_seqs_per_codon_groupi species in group 1 and $num_seqs_per_codon_groupj species in group 2\n";
#									}

                                  FIND_CONSERVED_CODON:
                                    for (
                                        my $gi_spp_index = 0 ;
                                        $gi_spp_index <
                                        $num_seqs_per_codon_groupi ;
                                        $gi_spp_index++
                                      )
                                    {

                                      FIND_CONSERVED_G2:
                                        for (
                                            my $gj_spp_index = 0 ;
                                            $gj_spp_index <
                                            $num_seqs_per_codon_groupj ;
                                            $gj_spp_index++
                                          )
                                        {   # for each gj sequence at this codon

                                            my $group_1_species =
                                              $group_1_spp[$gi_spp_index];
                                            my $group_2_species =
                                              $group_2_spp[$gj_spp_index];

                                            my $codon_gi_rep =
                                              $data_codonnum_spp_codon_hh{
                                                $codon_num}->{$group_1_species};
                                            my $codon_gj_rep =
                                              $data_codonnum_spp_codon_hh{
                                                $codon_num}->{$group_2_species};

#											if($codon_num == 2 && $partition_id == 1) {
#												print "\ncodon 2\ngroup_1_species=$group_1_species codon=$codon_gi_rep\ngroup_2_species=$group_2_species codon=$codon_gj_rep\n";
#											}

                                # grab first codon THAT DOESN'T CONTAIN N or GAP
                                            if (   $codon_gi_rep =~ 'N'
                                                || $codon_gi_rep =~ '-'
                                                || $codon_gi_rep eq '' )
                                            {
                                                if (   $codon_gj_rep =~ 'N'
                                                    || $codon_gj_rep =~ '-'
                                                    || $codon_gj_rep eq '' )
                                                {
#print "\ncodon $codon_num codons are codon_gi_rep $codon_gi_rep and codon_gj_rep $codon_gj_rep and neither work\n";
                                                    next FIND_CONSERVED_G2;
                                                }
                                                else {
                                                    $conserved_codon =
                                                      $codon_gj_rep;

#													if($codon_num == 2 && $partition_id == 1) {
#														print "\ncodon $codon_num codons are codon_gi_rep $codon_gi_rep and codon_gj_rep $codon_gj_rep and the last works\n";
#													}
                                                    last FIND_CONSERVED_CODON;
                                                }
                                            }
                                            else {
                                                $conserved_codon =
                                                  $codon_gi_rep;

#												if($codon_num == 2 && $partition_id == 1) {
#													print "\ncodon $codon_num codons are codon_gi_rep $codon_gi_rep and codon_gj_rep $codon_gj_rep and the first works\n";
#												}
                                                last FIND_CONSERVED_CODON;
                                            }
                                        }
                                    }

                               # Print a warning if no conserved codon was found
                                    if ( $conserved_codon eq '' ) {
                                        warn
"\n### WARNING: no codons found in between-group analysis at codon position $codon_num\.\n"
                                          . "### Most likely, the selected groups consist of all gaps at this codon, which is not necessarily a problem.\n"
                                          . "### Less likely, there may be multiple genes from a single species in family $family_id partition $partition_id.\n"
                                          . "### Inserting gaps and proceeding.\n\n";

                                        $conserved_codon = '---';
                                    }

                   #STORE THE CONSERVED CODON FOR LATER USE IN BOOTSTRAPPING
                   #					$codonIndex_conserved{$codon_index} = $conserved_codon;
                   #HERE

                                    my @codon_sites_1_arr =
                                      &get_number_of_sites( $conserved_codon,
                                        1 );
                                    my @codon_sites_2_arr =
                                      &get_number_of_sites( $conserved_codon,
                                        2 );
                                    my @codon_sites_3_arr =
                                      &get_number_of_sites( $conserved_codon,
                                        3 );

                                    my $codon_N_sites =
                                      ( $codon_sites_1_arr[0] +
                                          $codon_sites_2_arr[0] +
                                          $codon_sites_3_arr[0] );
                                    my $codon_S_sites =
                                      ( $codon_sites_1_arr[1] +
                                          $codon_sites_2_arr[1] +
                                          $codon_sites_3_arr[1] );

                         #					print CODON_FILE "conserved\t$conserved_codon\t".
                         #						"$codon_N_sites\t$codon_S_sites\t0\t0\n";

                                    $this_codon_output .=
                                      "conserved\t"
                                      . $codonum_group_numDefinedCodons_hh{
                                        $codon_num}->{group1}
                                      . "\t"
                                      . $codonum_group_numDefinedCodons_hh{
                                        $codon_num}->{group2}
                                      . "\t$conserved_codon\t"
                                      . "$codon_N_sites\t$codon_S_sites\t0\t0";

                                    # SLIDING WINDOW & BOOTSTRAP
                                    $codon_data_hh{$codon_num}->{N_sites} =
                                      $codon_N_sites;
                                    $codon_data_hh{$codon_num}->{S_sites} =
                                      $codon_S_sites;

                                   #					$product_N_sites_sum += $codon_N_sites;
                                   #					$product_S_sites_sum += $codon_S_sites;

#if($partition_id =~ /2K/) {
#	print "\nProduct $partition_id, conserved codon $conserved_codon N sites $codon_N_sites\n";
#}

                                }    # end not polymorphic

                                # ADD TO BUFFER
                                $output_buffer .= "$this_codon_output\n";
                                push( @codon_lines_arr, $this_codon_output );

                                # FLUSH BUFFER
                                if ( length($output_buffer) > 50000 ) {
                                    open( BETWEEN_CODON_FILE,
                                        ">>between\_group\_codon\_results\.txt"
                                    );
                                    print BETWEEN_CODON_FILE "$output_buffer";
                                    close BETWEEN_CODON_FILE;
                                    $output_buffer = '';
                                }

                            }    # end all codons in product

                            # FLUSH BUFFER
                            open( BETWEEN_CODON_FILE,
                                ">>between\_group\_codon\_results\.txt" );
                            print BETWEEN_CODON_FILE "$output_buffer";
                            close BETWEEN_CODON_FILE;
                            $output_buffer = '';

                            # Calculate overall values for the product
                            my $between_product_N_sites_sum = 0;
                            my $between_product_S_sites_sum = 0;
                            my $between_product_N_diffs_sum = 0;
                            my $between_product_S_diffs_sum = 0;

                            for (
                                my $line_index = 0 ;
                                $line_index < @codon_lines_arr ;
                                $line_index++
                              )
                            {
                                my $this_line = $codon_lines_arr[$line_index];
                                my @this_line_arr =
                                  split( /\t/, $this_line, -1 );
                                ##product / group_1 / group_2 / codon / variability / comparisons / N_sites / S_sites / N_diffs / S_diffs
#comparison1 / 80 / 1 /
#	(Apleuropneumoniae10D13039 Apleuropneumoniae10D13039 Aminor202 Aminor202 Asuccinogenes130Z Asuccinogenes130Z AactinomycetD11S1 AaphrophilusNJ8700 GanatisUMN179 GanatisUMN179 HaegyptiusATCC11116) /
#	(Hinfluenzae22121 Hinfluenzae22121 HinfluenzaeRdKW20 HparainfATCC33392 HparainfATCC33392 Hparasuis29755 Hparasuis29755 Hsomnus129PT MhaemolyticaA2BOV MhaemolyticaA2BOV Msuccinici Msuccinici PdagmatisATCC43325 PdagmatisATCC43325 PmultocidaPm70 PmultocidaPm70 EcoliK12MG1655 EcoliK12MG1655 EcoliK12MG1655  EcoliK12MG1655 EcoliK12MG1655) /
# codon_1  / conserved /  / 0 / 0 / 0 / 0

#analysis / family / partition / group_1 / group_2 / codon_num / variability / NUM_DEFINED_G1 / NUM_DEFINED_G2 / comparisons / N_sites / S_sites / N_diffs / S_diffs
                                my $poly = $this_line_arr[6];

                                # IF POLY
                                if ( $poly eq 'polymorphic' ) {
                                    $between_product_N_sites_sum +=
                                      $this_line_arr[10];
                                    $between_product_S_sites_sum +=
                                      $this_line_arr[11];
                                    $between_product_N_diffs_sum +=
                                      $this_line_arr[12];
                                    $between_product_S_diffs_sum +=
                                      $this_line_arr[13];

                                    # IF CONSERVED
                                }
                                elsif ( $poly eq 'conserved' ) {
                                    $between_product_N_sites_sum +=
                                      $this_line_arr[10];
                                    $between_product_S_sites_sum +=
                                      $this_line_arr[11];

                                    # PROBLEM
                                }
                                else {
                                    die
"\nPROBLEM: site neither conserved nor polymorphic!? DIE\n\n";
                                }

                            }    # END CODONS DATA

                            # Calculate overall values for the product
                            my $between_product_dN;
                            if ( $between_product_N_sites_sum > 0 ) {
                                $between_product_dN =
                                  $between_product_N_diffs_sum /
                                  $between_product_N_sites_sum;
                            }
                            else {
                                $between_product_dN = '*';
                            }

                            my $between_product_dS;
                            if ( $between_product_S_sites_sum > 0 ) {
                                $between_product_dS =
                                  $between_product_S_diffs_sum /
                                  $between_product_S_sites_sum;
                            }
                            else {
                                $between_product_dS = '*';
                            }

                            my $between_product_w;
                            if (   $between_product_dS > 0
                                && $between_product_dS ne '*' )
                            {
                                $between_product_w =
                                  $between_product_dN / $between_product_dS;
                            }
                            else {
                                $between_product_w = '*';
                            }

                            my $between_product_dN_minus_dS;
                            if (   $between_product_dN >= 0
                                && $between_product_dS >= 0 )
                            {
                                $between_product_dN_minus_dS =
                                  $between_product_dN - $between_product_dS;
                            }

                            $between_total_N_sites +=
                              $between_product_N_sites_sum;
                            $between_total_S_sites +=
                              $between_product_S_sites_sum;
                            $between_total_N_diffs +=
                              $between_product_N_diffs_sum;
                            $between_total_S_diffs +=
                              $between_product_S_diffs_sum;

        # Print product totals
        #my $out_line = "$i\t$j\t$partition_id\t$between_product_N_sites_sum\t".

#"analysis\tfamily\tpartition\tgroup_1\tgroup_2\tN_sites\tS_sites\tN_diffs\tS_diffs\t".
#	"dN\tdS\tdN_minus_dS\tdN_over_dS\n";

                            my $out_line =
"$analysis_name\t$family_id\t$partition_id\t$group_name_i\t$group_name_j\t$between_product_N_sites_sum\t"
                              . "$between_product_S_sites_sum\t$between_product_N_diffs_sum\t$between_product_S_diffs_sum\t"
                              . "$between_product_dN\t$between_product_dS\t$between_product_dN_minus_dS\t$between_product_w";

#							open(BETWEEN_OUTFILE,">>family$family_id\_between\_group\_results\.txt");
                            open( BETWEEN_OUTFILE,
                                ">>between\_group\_results\.txt" );

                            print BETWEEN_OUTFILE "$out_line\n";

                            #print "$out_line";

                            close BETWEEN_OUTFILE;

##							print "\nTotals:\nN_sites: $total_N_sites\nS_sites: $total_S_sites\n".
##								"N_diffs: $total_N_diffs\nS_diffs: $total_S_diffs\n\n";

                        }    # end group i vs. j

                        close BETWEEN_GROUP_CONFIG;

                    }
                    else {    # no between-group config file exists

                        warn
"\n### WARNING: between-group analysis requested but no file given. Skipping analysis.\n";

                    }

                }    # between-group flag was given

                # Print product totals
                my $product_piN;
                my $product_piS;
                my $product_piN_over_piS;

                if ( $product_N_sites_sum > 0 ) {
                    $product_piN = $product_N_diffs_sum / $product_N_sites_sum;
                }
                else {
                    $product_piN = '*';
                }

                if ( $product_S_sites_sum > 0 ) {
                    $product_piS = $product_S_diffs_sum / $product_S_sites_sum;
                }
                else {
                    $product_piS = '*';
                }

                my $product_piN_minus_piS = $product_piN - $product_piS;

                if ( $product_piS > 0 ) {
                    $product_piN_over_piS = $product_piN / $product_piS;
                }
                else {
                    $product_piN_over_piS = '*';
                }

                my $out_line =
                    "$family_id\t$partition_id\t$product_N_sites_sum\t"
                  . "$product_S_sites_sum\t$product_N_diffs_sum\t$product_S_diffs_sum\t"
                  . "$product_piN\t$product_piS\t$product_piN_minus_piS\t$product_piN_over_piS";

                $out_line .= "\n";

                #				chdir("$OID_USER_DIR\/data");
                #				open(BETWEEN_OUTFILE,">>family\_diversity\_results\.txt");
                #				print BETWEEN_OUTFILE "$out_line";
                #				#print "ACTUAL:\n$out_line\n\n";
                #				close BETWEEN_OUTFILE;
                #				chdir("$OID_USER_DIR\/data\/$family_id");

                open( WITHIN_GROUP_OUTFILE, ">>within_group_results.txt" );
                print WITHIN_GROUP_OUTFILE "$out_line";
                close WITHIN_GROUP_OUTFILE;

                #		open(VARIANT_FILE,">>within\_family\_variant\_results\.txt");

                # LOOP POLY CODONS; CODON NUMBER WILL BE INDEX + 1
                my %min_count2num_alleles;

                foreach
                  my $codon_index ( sort { $a <=> $b } keys %poly_site_data_hh )
                {
 # Here's what we did: $poly_site_data_hh{$codon_index}->{$codon_si} += $weight;
 # $weight = $comps_hh{$codon_si}->{$codon_sj};
 # $comps_hh{$codon_si}->{$codon_sj}+=1;
 # my $codon_si = $codonum_codon_aa[$codon_index]->[$si_seq_index];
 # push(@{$codonum_codon_aa[$array_index]},$codon);
 # my $num_seqs = scalar @{$codonum_codon_aa[0]};

                    my $codon_num = $codon_index + 1;

                    my $num_alleles =
                      scalar( keys %{ $poly_site_data_hh{$codon_index} } );

                    my $codons = '';

                    #			my $num_alleles = 0;
                    my $codon_counts = '';
                    my @codon_counts;
                    my $count_total = 0;

                    my $max_codon_count = 0;
                    my $consensus_codon = '';
                    my $ambiguous_flag  = 0;

                    foreach my $codon_triplet (
                        sort keys %{ $poly_site_data_hh{$codon_index} } )
                    {
                        $codons .= "$codon_triplet\,";

                        #				$num_alleles ++;
                        my $curr_codon_count;

                        if ( ( $num_seqs_defined[$codon_index] - 1 ) > 0 ) {
                            $curr_codon_count =
                              ( $poly_site_data_hh{$codon_index}
                                  ->{$codon_triplet} /
                                  ( $num_seqs_defined[$codon_index] - 1 ) );
                        }

#push(@codon_counts, ($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs_defined[$codon_index]-1)));
                        push( @codon_counts, $curr_codon_count );

#$count_total += ($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs_defined[$codon_index]-1));
                        $count_total += $curr_codon_count;
                        $codon_counts .= $curr_codon_count . "\,";

#				$codons .= "$codon_triplet\,";
#				$num_alleles ++;
#				my $curr_codon_count = ($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs-1));
#				$codon_counts .= $curr_codon_count . "\,";
#				push(@codon_counts,($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs-1)));
#				$count_total += ($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs-1));

                        # CONSENSUS CODON
                        if ( $curr_codon_count > $max_codon_count ) {
                            $max_codon_count = $curr_codon_count;
                            $consensus_codon = $codon_triplet;
                        }
                        elsif ( $curr_codon_count == $max_codon_count ) {
                            $ambiguous_flag = 1;
                        }
                    }

                    #my @freqs_arr;
                    my $codon_freqs = '';

                    foreach (@codon_counts) {
                        my $freq = ( $_ / $count_total );

                        #push(@freqs_arr, $freq);
                        my $freq_rounded = sprintf( "%.3f", $freq );
                        $codon_freqs .= "$freq_rounded\,";
                    }

                    # GOBBLE / TRIM ENDS OF EACH COMPONENT
                    chop($codons);
                    chop($codon_counts);
                    chop($codon_freqs);

                    # NUM SEQS MINOR ALLELE
                    my $min_codon_count =
                      $num_seqs_defined[$codon_index] - $max_codon_count;

                    #			my $min_codon_count = $num_seqs - $max_codon_count;

                    my $out_variant_line =
"$family_id\t$partition_id\t$codon_num\t$consensus_codon\t$num_alleles\t"
                      . "$codons\t$codon_counts\t$count_total\t$codon_freqs\t$max_codon_count\t$min_codon_count\t"
                      . $num_seqs_defined[$codon_index];

                    #			print VARIANT_FILE "$out_variant_line\n";

                    $min_count2num_alleles{$min_codon_count}->{$num_alleles}++;
                    $min_count2num_alleles_TOTALS{$min_codon_count}
                      ->{$num_alleles}++;

                }    # end loop poly codons within this product

                #		close VARIANT_FILE;

                $total_N_sites += $product_N_sites_sum;
                $total_S_sites += $product_S_sites_sum;
                $total_N_diffs += $product_N_diffs_sum;
                $total_S_diffs += $product_S_diffs_sum;

                foreach my $min_codon_count ( sort { $a <=> $b }
                    keys %min_count2num_alleles )
                {

                    foreach my $num_alleles ( sort { $a <=> $b }
                        keys %{ $min_count2num_alleles{$min_codon_count} } )
                    {

       #print "For codons with MAF of $min_codon_count, $num_alleles alleles: ".
       #	$min_count2num_alleles{$min_codon_count}->{$num_alleles} . " cases\n";

                    }
                }

                #close BETWEEN_OUTFILE;
                #					my $family_end_time = time;
                #					my $partition_end_time = time;

#					print FAMILY_LOG scalar(@sorted_aa_seqs_hash_keys) . "\t$partition_start_time\t$partition_end_time\n";
#					close FAMILY_LOG;
            }

            my $partition_end_time = time;

            print FAMILY_LOG scalar(@sorted_aa_seqs_hash_keys)
              . "\t$partition_start_time\t$partition_end_time\n";
            close FAMILY_LOG;

            print "\nDONE with partition.\n";

        }    # end all within-family partitions

    }
    else
    { # We've made sure that only directories with FAMILY.aligned files in them have been used
        print "No FAMILY.aligned file.\n";

        open( IGNOREFILE, ">>$OID_USER_DIR/data/$family_id/.snpgenie_ignore" );
        print IGNOREFILE "done";
        close IGNOREFILE;
    }
}

print "\nDONE with family.\n";

open( DONEFILE, ">>$OID_USER_DIR/data/$family_id/.snpgenie_done" );
print DONEFILE "done";
close DONEFILE;

#### Don't print the header if already present
###unless(-f "$OID_USER_DIR/data/family_diversity_results.txt") {
###	open(BETWEEN_OUTFILE,">>$OID_USER_DIR/data/family_diversity_results.txt");
###	print BETWEEN_OUTFILE "family\tpartition\tN_sites\tS_sites\tN_diffs\tS_diffs\tdN\tdS\tdN-dS\t".
###		"dN\/dS";
###
###	if($num_bootstraps > 1) {
###		print BETWEEN_OUTFILE "\tSE_dN-dS\tz_value\n";
###	} else {
###		print BETWEEN_OUTFILE "\n";
###	}
###	close BETWEEN_OUTFILE;
###}

###open(BETWEEN_OUTFILE,">>$OID_USER_DIR/data/family_diversity_results.txt");
###
###SWEEP: foreach my $family_id (@OID_USER_DIR_contents) {
###	if(-f "$OID_USER_DIR/data/$family_id/temp\_family\_diversity\_results\_line.txt") {
###		open(FAMILY_SUM, "$OID_USER_DIR/data/$family_id/temp\_family\_diversity\_results\_line.txt") or
###			die "Could not open file $OID_USER_DIR/data/$family_id/temp\_family\_diversity\_results\_line.txt";
###
###		while(<FAMILY_SUM>) {
###			#chomp;
###			print BETWEEN_OUTFILE "$_";
###		}
###
###		close FAMILY_SUM;
###		unlink "$OID_USER_DIR/data/$family_id/temp\_family\_diversity\_results\_line.txt";
###	}
###}
###close BETWEEN_OUTFILE;

##########################################################################################
## NEW NONCODING DIVERSITY FILE
#my @nc_A_count; # indices are alignment index (pos - 1)
#my @nc_C_count;
#my @nc_G_count;
#my @nc_T_count;
#
#my $aln_length = length($seqs_arr[0]);
#
## Initialize
#for(my $seq_pos = 0; $seq_pos < $aln_length; $seq_pos++) {
#	$nc_A_count[$seq_pos] = 0;
#	$nc_C_count[$seq_pos] = 0;
#	$nc_G_count[$seq_pos] = 0;
#	$nc_T_count[$seq_pos] = 0;
#}
#
#my $num_seqs = scalar(@seqs_arr);
#
## Store nucleotide counts
#for(my $seq_index = 0; $seq_index < $num_seqs; $seq_index++) {
#	my $this_seq = $seqs_arr[$seq_index];
#
#	for(my $seq_pos = 0; $seq_pos < $aln_length; $seq_pos++) {
#		my $this_nt = substr($this_seq,$seq_pos,1);
#		if($this_nt eq 'A') {
#			$nc_A_count[$seq_pos]++;
#		} elsif($this_nt eq 'C') {
#			$nc_C_count[$seq_pos]++;
#		} elsif($this_nt eq 'G') {
#			$nc_G_count[$seq_pos]++;
#		} elsif($this_nt eq 'T') {
#			$nc_T_count[$seq_pos]++;
#		}
#	}
#}
#
#my @pi_by_site;
#
## Calculate and store non-coding (just plain) pi
#open(BETWEEN_OUTFILE_GENE_DIV,">>within\_family\_site\_results\.txt");
#
#print BETWEEN_OUTFILE_GENE_DIV "file\tsite\tmaj_nt\tmaj_nt_count\t".
#	"pi\tcoverage\t".
#	"A\tC\tG\tT\n";
#
#for(my $seq_pos = 0; $seq_pos < $aln_length; $seq_pos++) {
#	my $A = $nc_A_count[$seq_pos];
#	my $C = $nc_C_count[$seq_pos];
#	my $G = $nc_G_count[$seq_pos];
#	my $T = $nc_T_count[$seq_pos];
#
#	my $maj_nt = 'A';
#	my $maj_nt_count = $A;
#
#	if($C > $maj_nt_count) {
#		$maj_nt = 'C';
#		$maj_nt_count = $C;
#	}
#
#	if($G > $maj_nt_count) {
#		$maj_nt = 'G';
#		$maj_nt_count = $G;
#	}
#
#	if($T > $maj_nt_count) {
#		$maj_nt = 'T';
#		$maj_nt_count = $T;
#	}
#
#	my $num_pw_diffs = ($A * $C) + ($A * $G) + ($A * $T) +
#		($C * $G) + ($C * $T) + ($G * $T);
#
#	my $num_pw_comps = ($num_seqs**2 - $num_seqs) / 2;
#
#	my $mean_num_pw_diffs = $num_pw_diffs / $num_pw_comps; # this is pi
#
#	$pi_by_site[$seq_pos] = $mean_num_pw_diffs;
#
#	my $curr_site = $seq_pos + 1;
#
#	print BETWEEN_OUTFILE_GENE_DIV "$nt_fasta_file_name\t$curr_site\t".
#		"$maj_nt\t$maj_nt_count\t".
#		"$mean_num_pw_diffs\t$num_seqs\t".
#		"$A\t$C\t$G\t$T\n";
#}
#close BETWEEN_OUTFILE_GENE_DIV;
#########################################################################################

# To track totals, no biological meaning for OrthologID
#print "\nTotals:\nN_sites: $total_N_sites\nS_sites: $total_S_sites\n".
#	"N_diffs: $total_N_diffs\nS_diffs: $total_S_diffs\n\n";

print "\n\n";

# Print a completion message to screen
&end_the_program;

#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################

#########################################################################################
sub median {
    my @values = sort { $a <=> $b } @_;
    my $length = @values;
    if ( $length % 2 ) {    # odd number of elements: return middle
        return $values[ int( $length / 2 ) ];
    }
    else {                  # even number of elements: return mean of middle two
        return (
            $values[ int( $length / 2 ) - 1 ] + $values[ int( $length / 2 ) ] )
          / 2;
    }
}

#########################################################################################
sub mean {
    my @values = @_;
    my $length = @values;
    my $sum;

    foreach (@values) {
        $sum += $_;
    }

    return ( $sum / $length );
}

#########################################################################################
sub standard_deviation {
    my @values         = @_;
    my $length         = @values;
    my $mean_of_values = &mean(@values);
    my $sum_squared_deviations;

    foreach (@values) {
        $sum_squared_deviations += ( $_ - $mean_of_values )**2;
    }

    my $variance = ($sum_squared_deviations) / ( $length - 1 );

    return ( sqrt($variance) );

}

#########################################################################################
sub align_codon2aa {

    #	my ($curr_aa_seq,$curr_nt_seq,$curr_seq_name) = @_;
    my ( $curr_aa_seq, $curr_nt_seq ) = @_;

# RECORD positions of GAPS in the amino acid, and also check for correct lengths
# Store gap indices for each seq in an array of arrays
    my @aa_gap_indices;

    my $curr_aa_seq_length          = length($curr_aa_seq);
    my $ungapped_curr_aa_seq_length = $curr_aa_seq_length;
    my $nt_seq_length               = length($curr_nt_seq);

    for ( my $pos_index = 0 ; $pos_index < $curr_aa_seq_length ; $pos_index++ )
    {
        my $curr_aa = substr( $curr_aa_seq, $pos_index, 1 );

        if ( $curr_aa =~ /-/ ) {    # It's a gap
            $ungapped_curr_aa_seq_length -= 1;
            push( @aa_gap_indices, $pos_index );
        }
    }

    if ( ( $ungapped_curr_aa_seq_length * 3 ) != $nt_seq_length ) {

# If nucleotide sequence contains a STOP codon not represented in the amino acid
# sequence or vice versa, modify the nucleotide sequence
        my $last_codon = substr( $curr_nt_seq, -3 );
        my $last_aa    = substr( $curr_aa_seq, -1 );
        if (
            (
                   $last_codon eq 'TAA'
                || $last_codon eq 'TAG'
                || $last_codon eq 'TGA'
            )
            && ( $last_aa ne '*' )
          )
        {
            # Remove the STOP codon
            substr( $curr_nt_seq, -3, 3, '' );
        }
        elsif (
            (
                   $last_codon ne 'TAA'
                && $last_codon ne 'TAG'
                && $last_codon ne 'TGA'
            )
            && ( $last_aa eq '*' )
          )
        {
            $curr_nt_seq .= 'TAA';
        }

        $nt_seq_length = length($curr_nt_seq);

        # If that didn't fix it
        if ( ( $ungapped_curr_aa_seq_length * 3 ) != $nt_seq_length ) {

#			die "\n\nDIE: A protein sequence does not contain " .
#				"the expected number of amino acids based on the length of nucleotide " .
#				"sequence. TERMINATED.\n\n";
#			print "\n\nWARNING: The protein sequence $curr_seq_name does not contain " .
#				"the expected number of amino acids based on the length of nucleotide " .
#				"sequence. Sequence populated with 'N'; excluded from dN/dS analysis.\n\n";

            return;

         #			my $seq_of_Ns;
         #
         #			for(my $counter=1; $counter<=$curr_aa_seq_length * 3; $counter++) {
         #				$seq_of_Ns .= "N";
         #			}
         #
         #			return $seq_of_Ns;
        }
    }

    # INSERT the appropriate gaps into the nucleotide sequences
    if ( scalar(@aa_gap_indices) >= 1 ) {

        foreach (@aa_gap_indices) {

            #print "$_ and ";
            substr( $curr_nt_seq, ( $_ * 3 ), 0, '---' );
        }

        #print "that's it.\n";

    #print "New nt sequence length of <THIS> is " . length($curr_nt_seq) . "\n";

    }

    # WRITE the new nucleotide sequence with amino acid alignment imposed
    my $curr_length = length($curr_nt_seq);

    if ( $curr_length < 3 ) {
        die "\n\n### WARNING: Must be at least 3 nucleotides.\n\n";
    }

    if ( ( $curr_length % 3 ) > 0 ) {
        warn
"\n\n### WARNING: Length is not a multiple of 3 (complete codons).\n\n";
    }

    #print "len is ".length($curr_nt_seq)."\n$curr_nt_seq";

    #print "\nAligned Nucleotide Sequence:\n$curr_seq\n\n";

    return $curr_nt_seq;

}

#########################################################################################
# Print a FASTA file to working directory with the given header and sequence
sub print_fasta_file {
    my ( $name, $header, $sequence ) = @_;
    open( OUT, ">>$name" );
    my $seq_length = length($sequence);

    print OUT "$header\n";

    for ( my $j = 0 ; $j < $seq_length ; $j += 60 ) {
        if ( $j > ( $seq_length - 60 ) ) {
            my $line = substr( $sequence, $j );
            print OUT "$line\n";
            last;
        }
        else {
            my $line = substr( $sequence, $j, 60 );
            print OUT "$line\n";
        }
    }

    close OUT;
}

#########################################################################################
# Get the number of synonymous and nonsynonymous sites in a particular codon position
sub get_number_of_sites {
    my ( $codon, $position ) = @_;

    $codon = uc($codon);
    $codon =~ tr/U/T/;

    my $codon_N_position = "$codon" . "_N_" . "$position";
    my $codon_S_position = "$codon" . "_S_" . "$position";
    my $num_N;
    my $num_S;
    my @num_sites;

# Numbers of nonsynonymous and synonymous sites for all three sites of every codon.
# Per Nei-Gojobori, STOP codons have 0 sites. Results in arrays codon name in the format:
# SITE1N SITE1S SITE2N SITE2S SITE3N SITE3S
    my %codon_type_position_hash = (
        "AAA_N_1" => 1,
        "AAA_S_1" => 0,
        "AAA_N_2" => 1,
        "AAA_S_2" => 0,
        "AAA_N_3" => 2 / 3,
        "AAA_S_3" => 1 / 3,
        "AAC_N_1" => 1,
        "AAC_S_1" => 0,
        "AAC_N_2" => 1,
        "AAC_S_2" => 0,
        "AAC_N_3" => 2 / 3,
        "AAC_S_3" => 1 / 3,
        "AAG_N_1" => 1,
        "AAG_S_1" => 0,
        "AAG_N_2" => 1,
        "AAG_S_2" => 0,
        "AAG_N_3" => 2 / 3,
        "AAG_S_3" => 1 / 3,
        "AAT_N_1" => 1,
        "AAT_S_1" => 0,
        "AAT_N_2" => 1,
        "AAT_S_2" => 0,
        "AAT_N_3" => 2 / 3,
        "AAT_S_3" => 1 / 3,
        "ACA_N_1" => 1,
        "ACA_S_1" => 0,
        "ACA_N_2" => 1,
        "ACA_S_2" => 0,
        "ACA_N_3" => 0,
        "ACA_S_3" => 1,
        "ACC_N_1" => 1,
        "ACC_S_1" => 0,
        "ACC_N_2" => 1,
        "ACC_S_2" => 0,
        "ACC_N_3" => 0,
        "ACC_S_3" => 1,
        "ACG_N_1" => 1,
        "ACG_S_1" => 0,
        "ACG_N_2" => 1,
        "ACG_S_2" => 0,
        "ACG_N_3" => 0,
        "ACG_S_3" => 1,
        "ACT_N_1" => 1,
        "ACT_S_1" => 0,
        "ACT_N_2" => 1,
        "ACT_S_2" => 0,
        "ACT_N_3" => 0,
        "ACT_S_3" => 1,
        "AGA_N_1" => 1 / 2,
        "AGA_S_1" => 1 / 2,
        "AGA_N_2" => 1,
        "AGA_S_2" => 0,
        "AGA_N_3" => 2 / 3,
        "AGA_S_3" => 1 / 3,
        "AGC_N_1" => 1,
        "AGC_S_1" => 0,
        "AGC_N_2" => 1,
        "AGC_S_2" => 0,
        "AGC_N_3" => 2 / 3,
        "AGC_S_3" => 1 / 3,
        "AGG_N_1" => 2 / 3,
        "AGG_S_1" => 1 / 3,
        "AGG_N_2" => 1,
        "AGG_S_2" => 0,
        "AGG_N_3" => 2 / 3,
        "AGG_S_3" => 1 / 3,
        "AGT_N_1" => 1,
        "AGT_S_1" => 0,
        "AGT_N_2" => 1,
        "AGT_S_2" => 0,
        "AGT_N_3" => 2 / 3,
        "AGT_S_3" => 1 / 3,
        "ATA_N_1" => 1,
        "ATA_S_1" => 0,
        "ATA_N_2" => 1,
        "ATA_S_2" => 0,
        "ATA_N_3" => 1 / 3,
        "ATA_S_3" => 2 / 3,
        "ATC_N_1" => 1,
        "ATC_S_1" => 0,
        "ATC_N_2" => 1,
        "ATC_S_2" => 0,
        "ATC_N_3" => 1 / 3,
        "ATC_S_3" => 2 / 3,
        "ATG_N_1" => 1,
        "ATG_S_1" => 0,
        "ATG_N_2" => 1,
        "ATG_S_2" => 0,
        "ATG_N_3" => 1,
        "ATG_S_3" => 0,
        "ATT_N_1" => 1,
        "ATT_S_1" => 0,
        "ATT_N_2" => 1,
        "ATT_S_2" => 0,
        "ATT_N_3" => 1 / 3,
        "ATT_S_3" => 2 / 3,

        "CAA_N_1" => 1,
        "CAA_S_1" => 0,
        "CAA_N_2" => 1,
        "CAA_S_2" => 0,
        "CAA_N_3" => 2 / 3,
        "CAA_S_3" => 1 / 3,
        "CAC_N_1" => 1,
        "CAC_S_1" => 0,
        "CAC_N_2" => 1,
        "CAC_S_2" => 0,
        "CAC_N_3" => 2 / 3,
        "CAC_S_3" => 1 / 3,
        "CAG_N_1" => 1,
        "CAG_S_1" => 0,
        "CAG_N_2" => 1,
        "CAG_S_2" => 0,
        "CAG_N_3" => 2 / 3,
        "CAG_S_3" => 1 / 3,
        "CAT_N_1" => 1,
        "CAT_S_1" => 0,
        "CAT_N_2" => 1,
        "CAT_S_2" => 0,
        "CAT_N_3" => 2 / 3,
        "CAT_S_3" => 1 / 3,
        "CCA_N_1" => 1,
        "CCA_S_1" => 0,
        "CCA_N_2" => 1,
        "CCA_S_2" => 0,
        "CCA_N_3" => 0,
        "CCA_S_3" => 1,
        "CCC_N_1" => 1,
        "CCC_S_1" => 0,
        "CCC_N_2" => 1,
        "CCC_S_2" => 0,
        "CCC_N_3" => 0,
        "CCC_S_3" => 1,
        "CCG_N_1" => 1,
        "CCG_S_1" => 0,
        "CCG_N_2" => 1,
        "CCG_S_2" => 0,
        "CCG_N_3" => 0,
        "CCG_S_3" => 1,
        "CCT_N_1" => 1,
        "CCT_S_1" => 0,
        "CCT_N_2" => 1,
        "CCT_S_2" => 0,
        "CCT_N_3" => 0,
        "CCT_S_3" => 1,
        "CGA_N_1" => 1 / 2,
        "CGA_S_1" => 1 / 2,
        "CGA_N_2" => 1,
        "CGA_S_2" => 0,
        "CGA_N_3" => 0,
        "CGA_S_3" => 1,
        "CGC_N_1" => 1,
        "CGC_S_1" => 0,
        "CGC_N_2" => 1,
        "CGC_S_2" => 0,
        "CGC_N_3" => 0,
        "CGC_S_3" => 1,
        "CGG_N_1" => 2 / 3,
        "CGG_S_1" => 1 / 3,
        "CGG_N_2" => 1,
        "CGG_S_2" => 0,
        "CGG_N_3" => 0,
        "CGG_S_3" => 1,
        "CGT_N_1" => 1,
        "CGT_S_1" => 0,
        "CGT_N_2" => 1,
        "CGT_S_2" => 0,
        "CGT_N_3" => 0,
        "CGT_S_3" => 1,
        "CTA_N_1" => 2 / 3,
        "CTA_S_1" => 1 / 3,
        "CTA_N_2" => 1,
        "CTA_S_2" => 0,
        "CTA_N_3" => 0,
        "CTA_S_3" => 1,
        "CTC_N_1" => 1,
        "CTC_S_1" => 0,
        "CTC_N_2" => 1,
        "CTC_S_2" => 0,
        "CTC_N_3" => 0,
        "CTC_S_3" => 1,
        "CTG_N_1" => 2 / 3,
        "CTG_S_1" => 1 / 3,
        "CTG_N_2" => 1,
        "CTG_S_2" => 0,
        "CTG_N_3" => 0,
        "CTG_S_3" => 1,
        "CTT_N_1" => 1,
        "CTT_S_1" => 0,
        "CTT_N_2" => 1,
        "CTT_S_2" => 0,
        "CTT_N_3" => 0,
        "CTT_S_3" => 1,

        "GAA_N_1" => 1,
        "GAA_S_1" => 0,
        "GAA_N_2" => 1,
        "GAA_S_2" => 0,
        "GAA_N_3" => 2 / 3,
        "GAA_S_3" => 1 / 3,
        "GAC_N_1" => 1,
        "GAC_S_1" => 0,
        "GAC_N_2" => 1,
        "GAC_S_2" => 0,
        "GAC_N_3" => 2 / 3,
        "GAC_S_3" => 1 / 3,
        "GAG_N_1" => 1,
        "GAG_S_1" => 0,
        "GAG_N_2" => 1,
        "GAG_S_2" => 0,
        "GAG_N_3" => 2 / 3,
        "GAG_S_3" => 1 / 3,
        "GAT_N_1" => 1,
        "GAT_S_1" => 0,
        "GAT_N_2" => 1,
        "GAT_S_2" => 0,
        "GAT_N_3" => 2 / 3,
        "GAT_S_3" => 1 / 3,
        "GCA_N_1" => 1,
        "GCA_S_1" => 0,
        "GCA_N_2" => 1,
        "GCA_S_2" => 0,
        "GCA_N_3" => 0,
        "GCA_S_3" => 1,
        "GCC_N_1" => 1,
        "GCC_S_1" => 0,
        "GCC_N_2" => 1,
        "GCC_S_2" => 0,
        "GCC_N_3" => 0,
        "GCC_S_3" => 1,
        "GCG_N_1" => 1,
        "GCG_S_1" => 0,
        "GCG_N_2" => 1,
        "GCG_S_2" => 0,
        "GCG_N_3" => 0,
        "GCG_S_3" => 1,
        "GCT_N_1" => 1,
        "GCT_S_1" => 0,
        "GCT_N_2" => 1,
        "GCT_S_2" => 0,
        "GCT_N_3" => 0,
        "GCT_S_3" => 1,
        "GGA_N_1" => 1,
        "GGA_S_1" => 0,
        "GGA_N_2" => 1,
        "GGA_S_2" => 0,
        "GGA_N_3" => 0,
        "GGA_S_3" => 1,
        "GGC_N_1" => 1,
        "GGC_S_1" => 0,
        "GGC_N_2" => 1,
        "GGC_S_2" => 0,
        "GGC_N_3" => 0,
        "GGC_S_3" => 1,
        "GGG_N_1" => 1,
        "GGG_S_1" => 0,
        "GGG_N_2" => 1,
        "GGG_S_2" => 0,
        "GGG_N_3" => 0,
        "GGG_S_3" => 1,
        "GGT_N_1" => 1,
        "GGT_S_1" => 0,
        "GGT_N_2" => 1,
        "GGT_S_2" => 0,
        "GGT_N_3" => 0,
        "GGT_S_3" => 1,
        "GTA_N_1" => 1,
        "GTA_S_1" => 0,
        "GTA_N_2" => 1,
        "GTA_S_2" => 0,
        "GTA_N_3" => 0,
        "GTA_S_3" => 1,
        "GTC_N_1" => 1,
        "GTC_S_1" => 0,
        "GTC_N_2" => 1,
        "GTC_S_2" => 0,
        "GTC_N_3" => 0,
        "GTC_S_3" => 1,
        "GTG_N_1" => 1,
        "GTG_S_1" => 0,
        "GTG_N_2" => 1,
        "GTG_S_2" => 0,
        "GTG_N_3" => 0,
        "GTG_S_3" => 1,
        "GTT_N_1" => 1,
        "GTT_S_1" => 0,
        "GTT_N_2" => 1,
        "GTT_S_2" => 0,
        "GTT_N_3" => 0,
        "GTT_S_3" => 1,

        "TAA_N_1" => 0,
        "TAA_S_1" => 0,
        "TAA_N_2" => 0,
        "TAA_S_2" => 0,
        "TAA_N_3" => 0,
        "TAA_S_3" => 0,       # STOP
        "TAC_N_1" => 1,
        "TAC_S_1" => 0,
        "TAC_N_2" => 1,
        "TAC_S_2" => 0,
        "TAC_N_3" => 0,
        "TAC_S_3" => 1,
        "TAG_N_1" => 0,
        "TAG_S_1" => 0,
        "TAG_N_2" => 0,
        "TAG_S_2" => 0,
        "TAG_N_3" => 0,
        "TAG_S_3" => 0,       # STOP
        "TAT_N_1" => 1,
        "TAT_S_1" => 0,
        "TAT_N_2" => 1,
        "TAT_S_2" => 0,
        "TAT_N_3" => 0,
        "TAT_S_3" => 1,
        "TCA_N_1" => 1,
        "TCA_S_1" => 0,
        "TCA_N_2" => 1,
        "TCA_S_2" => 0,
        "TCA_N_3" => 0,
        "TCA_S_3" => 1,
        "TCC_N_1" => 1,
        "TCC_S_1" => 0,
        "TCC_N_2" => 1,
        "TCC_S_2" => 0,
        "TCC_N_3" => 0,
        "TCC_S_3" => 1,
        "TCG_N_1" => 1,
        "TCG_S_1" => 0,
        "TCG_N_2" => 1,
        "TCG_S_2" => 0,
        "TCG_N_3" => 0,
        "TCG_S_3" => 1,
        "TCT_N_1" => 1,
        "TCT_S_1" => 0,
        "TCT_N_2" => 1,
        "TCT_S_2" => 0,
        "TCT_N_3" => 0,
        "TCT_S_3" => 1,
        "TGA_N_1" => 0,
        "TGA_S_1" => 0,
        "TGA_N_2" => 0,
        "TGA_S_2" => 0,
        "TGA_N_3" => 0,
        "TGA_S_3" => 0,       # STOP
        "TGC_N_1" => 1,
        "TGC_S_1" => 0,
        "TGC_N_2" => 1,
        "TGC_S_2" => 0,
        "TGC_N_3" => 1 / 2,
        "TGC_S_3" => 1 / 2,
        "TGG_N_1" => 1,
        "TGG_S_1" => 0,
        "TGG_N_2" => 1,
        "TGG_S_2" => 0,
        "TGG_N_3" => 1,
        "TGG_S_3" => 0,
        "TGT_N_1" => 1,
        "TGT_S_1" => 0,
        "TGT_N_2" => 1,
        "TGT_S_2" => 0,
        "TGT_N_3" => 1 / 2,
        "TGT_S_3" => 1 / 2,
        "TTA_N_1" => 2 / 3,
        "TTA_S_1" => 1 / 3,
        "TTA_N_2" => 1,
        "TTA_S_2" => 0,
        "TTA_N_3" => 2 / 3,
        "TTA_S_3" => 1 / 3,
        "TTC_N_1" => 1,
        "TTC_S_1" => 0,
        "TTC_N_2" => 1,
        "TTC_S_2" => 0,
        "TTC_N_3" => 2 / 3,
        "TTC_S_3" => 1 / 3,
        "TTG_N_1" => 2 / 3,
        "TTG_S_1" => 1 / 3,
        "TTG_N_2" => 1,
        "TTG_S_2" => 0,
        "TTG_N_3" => 2 / 3,
        "TTG_S_3" => 1 / 3,
        "TTT_N_1" => 1,
        "TTT_S_1" => 0,
        "TTT_N_2" => 1,
        "TTT_S_2" => 0,
        "TTT_N_3" => 2 / 3,
        "TTT_S_3" => 1 / 3,
    );

    $num_N = $codon_type_position_hash{$codon_N_position};
    $num_S = $codon_type_position_hash{$codon_S_position};

    @num_sites = ( $num_N, $num_S );

    return @num_sites;
}

#########################################################################################
# Returns the number of differences between two codons, averaged over all minimal-evolution paths
# Per Nei-Gojobori, we return 0 differences for STOP codons, because they have 0 sites.
sub return_avg_diffs {
    my ( $codon1, $codon2 ) = @_;

    my %all_diffs_hh = (
        'AAA' => {
            'AAC' => { 'N' => 1,                'S' => 0 },
            'AAG' => { 'N' => 0,                'S' => 1 },
            'AAT' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 1,                'S' => 0 },
            'ACC' => { 'N' => 1.5,              'S' => 0.5 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1.5,              'S' => 0.5 },
            'AGA' => { 'N' => 1,                'S' => 0 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 1,                'S' => 1 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 1,                'S' => 0 },
            'ATC' => { 'N' => 1.5,              'S' => 0.5 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1.5,              'S' => 0.5 },
            'CAA' => { 'N' => 1,                'S' => 0 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 1,                'S' => 1 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2,                'S' => 0 },
            'CCC' => { 'N' => 2.5,              'S' => 0.5 },
            'CCG' => { 'N' => 2,                'S' => 1 },
            'CCT' => { 'N' => 2.5,              'S' => 0.5 },
            'CGA' => { 'N' => 1.5,              'S' => 0.5 },
            'CGC' => { 'N' => 2.5,              'S' => 0.5 },
            'CGG' => { 'N' => 1.5,              'S' => 1.5 },
            'CGT' => { 'N' => 2.5,              'S' => 0.5 },
            'CTA' => { 'N' => 2,                'S' => 0 },
            'CTC' => { 'N' => 2.5,              'S' => 0.5 },
            'CTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTT' => { 'N' => 2.5,              'S' => 0.5 },
            'GAA' => { 'N' => 1,                'S' => 0 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 1,                'S' => 1 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2,                'S' => 0 },
            'GCC' => { 'N' => 2.5,              'S' => 0.5 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2.5,              'S' => 0.5 },
            'GGA' => { 'N' => 2,                'S' => 0 },
            'GGC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GTA' => { 'N' => 2,                'S' => 0 },
            'GTC' => { 'N' => 2.5,              'S' => 0.5 },
            'GTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTT' => { 'N' => 2.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 0 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 3,                'S' => 0 },
            'TGG' => { 'N' => 2,                'S' => 1 },
            'TGT' => { 'N' => 3,                'S' => 0 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 2.75,             'S' => 0.25 },
            'TTG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTT' => { 'N' => 2.75,             'S' => 0.25 },
        },
        'AAC' => {
            'AAA' => { 'N' => 1,                'S' => 0 },
            'AAG' => { 'N' => 1,                'S' => 0 },
            'AAT' => { 'N' => 0,                'S' => 1 },
            'ACA' => { 'N' => 1.5,              'S' => 0.5 },
            'ACC' => { 'N' => 1,                'S' => 0 },
            'ACG' => { 'N' => 1.5,              'S' => 0.5 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 1,                'S' => 0 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 1,                'S' => 1 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1,                'S' => 0 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 1,                'S' => 0 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 1,                'S' => 1 },
            'CCA' => { 'N' => 2.5,              'S' => 0.5 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2.5,              'S' => 0.5 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2.5,              'S' => 0.5 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 1,                'S' => 0 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 1,                'S' => 1 },
            'GCA' => { 'N' => 2.5,              'S' => 0.5 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2.5,              'S' => 0.5 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2.5,              'S' => 0.5 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 1,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 1,                'S' => 1 },
            'TCA' => { 'N' => 2.25,             'S' => 0.75 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2.25,             'S' => 0.75 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 3,                'S' => 0 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 2.75,             'S' => 0.25 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 3,                'S' => 0 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'AAG' => {
            'AAA' => { 'N' => 0,                'S' => 1 },
            'AAC' => { 'N' => 1,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1.5,              'S' => 0.5 },
            'ACG' => { 'N' => 1,                'S' => 0 },
            'ACT' => { 'N' => 1.5,              'S' => 0.5 },
            'AGA' => { 'N' => 1,                'S' => 1 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 1,                'S' => 0 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 1,                'S' => 1 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 1,                'S' => 0 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2.5,              'S' => 0.5 },
            'CCG' => { 'N' => 2,                'S' => 0 },
            'CCT' => { 'N' => 2.5,              'S' => 0.5 },
            'CGA' => { 'N' => 1.5,              'S' => 1.5 },
            'CGC' => { 'N' => 2.5,              'S' => 0.5 },
            'CGG' => { 'N' => 1.5,              'S' => 0.5 },
            'CGT' => { 'N' => 2.5,              'S' => 0.5 },
            'CTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CTG' => { 'N' => 2,                'S' => 0 },
            'CTT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GAA' => { 'N' => 1,                'S' => 1 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 1,                'S' => 0 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2.5,              'S' => 0.5 },
            'GCG' => { 'N' => 2,                'S' => 0 },
            'GCT' => { 'N' => 2.5,              'S' => 0.5 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GGG' => { 'N' => 2,                'S' => 0 },
            'GGT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GTG' => { 'N' => 2,                'S' => 0 },
            'GTT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 0 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 3,                'S' => 0 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 3,                'S' => 0 },
            'TTA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTC' => { 'N' => 3,                'S' => 0 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 3,                'S' => 0 },
        },
        'AAT' => {
            'AAA' => { 'N' => 1,                'S' => 0 },
            'AAC' => { 'N' => 0,                'S' => 1 },
            'AAG' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 1.5,              'S' => 0.5 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1.5,              'S' => 0.5 },
            'ACT' => { 'N' => 1,                'S' => 0 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 1,                'S' => 1 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 1,                'S' => 0 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 1,                'S' => 1 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 2.5,              'S' => 0.5 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2.5,              'S' => 0.5 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 2.5,              'S' => 0.5 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 1,                'S' => 1 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 2.5,              'S' => 0.5 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2.5,              'S' => 0.5 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 2.5,              'S' => 0.5 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 1,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 1,                'S' => 0 },
            'TCA' => { 'N' => 2.25,             'S' => 0.75 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2.25,             'S' => 0.75 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 3,                'S' => 0 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 2.75,             'S' => 0.25 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 3,                'S' => 0 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'ACA' => {
            'AAA' => { 'N' => 1,                'S' => 0 },
            'AAC' => { 'N' => 1.5,              'S' => 0.5 },
            'AAG' => { 'N' => 1,                'S' => 1 },
            'AAT' => { 'N' => 1.5,              'S' => 0.5 },
            'ACC' => { 'N' => 0,                'S' => 1 },
            'ACG' => { 'N' => 0,                'S' => 1 },
            'ACT' => { 'N' => 0,                'S' => 1 },
            'AGA' => { 'N' => 1,                'S' => 0 },
            'AGC' => { 'N' => 1.5,              'S' => 0.5 },
            'AGG' => { 'N' => 1,                'S' => 1 },
            'AGT' => { 'N' => 1.5,              'S' => 0.5 },
            'ATA' => { 'N' => 1,                'S' => 0 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2,                'S' => 1 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1,                'S' => 0 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 1.5,              'S' => 0.5 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 1.5,              'S' => 1.5 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 0 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2,                'S' => 1 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 1,                'S' => 0 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 0 },
            'GGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTA' => { 'N' => 2,                'S' => 0 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.25,             'S' => 0.75 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.25,             'S' => 0.75 },
            'TCA' => { 'N' => 1,                'S' => 0 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2.25,             'S' => 0.75 },
            'TGG' => { 'N' => 2,                'S' => 1 },
            'TGT' => { 'N' => 2.25,             'S' => 0.75 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TTT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
        },
        'ACC' => {
            'AAA' => { 'N' => 1.5,              'S' => 0.5 },
            'AAC' => { 'N' => 1,                'S' => 0 },
            'AAG' => { 'N' => 1.5,              'S' => 0.5 },
            'AAT' => { 'N' => 1,                'S' => 1 },
            'ACA' => { 'N' => 0,                'S' => 1 },
            'ACG' => { 'N' => 0,                'S' => 1 },
            'ACT' => { 'N' => 0,                'S' => 1 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 1,                'S' => 0 },
            'AGG' => { 'N' => 1.5,              'S' => 0.5 },
            'AGT' => { 'N' => 1,                'S' => 1 },
            'ATA' => { 'N' => 1,                'S' => 1 },
            'ATC' => { 'N' => 1,                'S' => 0 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 0 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 0 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 0 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 2.5,              'S' => 0.5 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 2.5,              'S' => 0.5 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'ACG' => {
            'AAA' => { 'N' => 1,                'S' => 1 },
            'AAC' => { 'N' => 1.5,              'S' => 0.5 },
            'AAG' => { 'N' => 1,                'S' => 0 },
            'AAT' => { 'N' => 1.5,              'S' => 0.5 },
            'ACA' => { 'N' => 0,                'S' => 1 },
            'ACC' => { 'N' => 0,                'S' => 1 },
            'ACT' => { 'N' => 0,                'S' => 1 },
            'AGA' => { 'N' => 1,                'S' => 1 },
            'AGC' => { 'N' => 1.5,              'S' => 0.5 },
            'AGG' => { 'N' => 1,                'S' => 0 },
            'AGT' => { 'N' => 1.5,              'S' => 0.5 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1.5,              'S' => 0.5 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 1.5,              'S' => 0.5 },
            'CAA' => { 'N' => 2,                'S' => 1 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 0 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 1.5,              'S' => 1.5 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 1.5,              'S' => 0.5 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTG' => { 'N' => 2,                'S' => 0 },
            'CTT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GAA' => { 'N' => 2,                'S' => 1 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 0 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGG' => { 'N' => 2,                'S' => 0 },
            'GGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTG' => { 'N' => 2,                'S' => 0 },
            'GTT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.25,             'S' => 0.75 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.25,             'S' => 0.75 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 0 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2.5,              'S' => 0.5 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 2.5,              'S' => 0.5 },
            'TTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TTC' => { 'N' => 2.5,              'S' => 0.5 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 2.5,              'S' => 0.5 },
        },
        'ACT' => {
            'AAA' => { 'N' => 1.5,              'S' => 0.5 },
            'AAC' => { 'N' => 1,                'S' => 1 },
            'AAG' => { 'N' => 1.5,              'S' => 0.5 },
            'AAT' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 0,                'S' => 1 },
            'ACC' => { 'N' => 0,                'S' => 1 },
            'ACG' => { 'N' => 0,                'S' => 1 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 1,                'S' => 1 },
            'AGG' => { 'N' => 1.5,              'S' => 0.5 },
            'AGT' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 1,                'S' => 1 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 0 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 0 },
            'CGA' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 0 },
            'GGA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2.5,              'S' => 0.5 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 2.5,              'S' => 0.5 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'AGA' => {
            'AAA' => { 'N' => 1,                'S' => 0 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 1,                'S' => 1 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 1,                'S' => 0 },
            'ACC' => { 'N' => 1.5,              'S' => 0.5 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 1,                'S' => 0 },
            'AGG' => { 'N' => 0,                'S' => 1 },
            'AGT' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 1,                'S' => 0 },
            'ATC' => { 'N' => 1.5,              'S' => 0.5 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1.5,              'S' => 0.5 },
            'CAA' => { 'N' => 1.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 1.5,              'S' => 1.5 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1.5,              'S' => 0.5 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 1.5,              'S' => 1.5 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 0,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 0,                'S' => 2 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 1.5,              'S' => 0.5 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 1.66666666666667, 'S' => 1.33333333333333 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAG' => { 'N' => 2,                'S' => 1 },
            'GAT' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GCA' => { 'N' => 2,                'S' => 0 },
            'GCC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGA' => { 'N' => 1,                'S' => 0 },
            'GGC' => { 'N' => 1.5,              'S' => 0.5 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1.5,              'S' => 0.5 },
            'GTA' => { 'N' => 2,                'S' => 0 },
            'GTC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 3,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 3,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 0 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 1,                'S' => 1 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 2.75,             'S' => 0.25 },
            'TTG' => { 'N' => 2.25,             'S' => 0.75 },
            'TTT' => { 'N' => 2.75,             'S' => 0.25 },
        },
        'AGC' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 1,                'S' => 0 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 1 },
            'ACA' => { 'N' => 1.5,              'S' => 0.5 },
            'ACC' => { 'N' => 1,                'S' => 0 },
            'ACG' => { 'N' => 1.5,              'S' => 0.5 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 1,                'S' => 0 },
            'AGG' => { 'N' => 1,                'S' => 0 },
            'AGT' => { 'N' => 0,                'S' => 1 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1,                'S' => 0 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 0 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 1.5,              'S' => 0.5 },
            'GGC' => { 'N' => 1,                'S' => 0 },
            'GGG' => { 'N' => 1.5,              'S' => 0.5 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2.5,              'S' => 0.5 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 2.25,             'S' => 0.75 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2.5,              'S' => 0.5 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 1,                'S' => 0 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 2.75,             'S' => 0.25 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 3,                'S' => 0 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'AGG' => {
            'AAA' => { 'N' => 1,                'S' => 1 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 1,                'S' => 0 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1.5,              'S' => 0.5 },
            'ACG' => { 'N' => 1,                'S' => 0 },
            'ACT' => { 'N' => 1.5,              'S' => 0.5 },
            'AGA' => { 'N' => 0,                'S' => 1 },
            'AGC' => { 'N' => 1,                'S' => 0 },
            'AGT' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 1.5,              'S' => 1.5 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 1.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1.5,              'S' => 1.5 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 1.5,              'S' => 0.5 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 0,                'S' => 2 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 0,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 1.66666666666667, 'S' => 1.33333333333333 },
            'CTC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTG' => { 'N' => 1.5,              'S' => 0.5 },
            'CTT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GAA' => { 'N' => 2,                'S' => 1 },
            'GAC' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCG' => { 'N' => 2,                'S' => 0 },
            'GCT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1.5,              'S' => 0.5 },
            'GGG' => { 'N' => 1,                'S' => 0 },
            'GGT' => { 'N' => 1.5,              'S' => 0.5 },
            'GTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTC' => { 'N' => 2.5,              'S' => 0.5 },
            'GTG' => { 'N' => 2,                'S' => 0 },
            'GTT' => { 'N' => 2.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 3,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 3,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 0 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 1,                'S' => 0 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 2.25,             'S' => 0.75 },
            'TTC' => { 'N' => 3,                'S' => 0 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 3,                'S' => 0 },
        },
        'AGT' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 1,                'S' => 1 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 1.5,              'S' => 0.5 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1.5,              'S' => 0.5 },
            'ACT' => { 'N' => 1,                'S' => 0 },
            'AGA' => { 'N' => 1,                'S' => 0 },
            'AGC' => { 'N' => 0,                'S' => 1 },
            'AGG' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 1,                'S' => 0 },
            'CAA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 0 },
            'CTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 1.5,              'S' => 0.5 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1.5,              'S' => 0.5 },
            'GGT' => { 'N' => 1,                'S' => 0 },
            'GTA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2.5,              'S' => 0.5 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2.25,             'S' => 0.75 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2.5,              'S' => 0.5 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 1,                'S' => 0 },
            'TTA' => { 'N' => 2.75,             'S' => 0.25 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 3,                'S' => 0 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'ATA' => {
            'AAA' => { 'N' => 1,                'S' => 0 },
            'AAC' => { 'N' => 1.5,              'S' => 0.5 },
            'AAG' => { 'N' => 1.5,              'S' => 0.5 },
            'AAT' => { 'N' => 1.5,              'S' => 0.5 },
            'ACA' => { 'N' => 1,                'S' => 0 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1.5,              'S' => 0.5 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 1,                'S' => 0 },
            'AGC' => { 'N' => 1.5,              'S' => 0.5 },
            'AGG' => { 'N' => 1.5,              'S' => 0.5 },
            'AGT' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 0,                'S' => 1 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 0,                'S' => 1 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 2,                'S' => 0 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1.5,              'S' => 0.5 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 0 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1.5,              'S' => 0.5 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 2,                'S' => 0 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 0 },
            'GGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTA' => { 'N' => 1,                'S' => 0 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1.5,              'S' => 0.5 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2.5,              'S' => 0.5 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2.5,              'S' => 0.5 },
            'TCA' => { 'N' => 2,                'S' => 0 },
            'TCC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TCT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2.5,              'S' => 0.5 },
            'TGG' => { 'N' => 2.5,              'S' => 0.5 },
            'TGT' => { 'N' => 2.5,              'S' => 0.5 },
            'TTA' => { 'N' => 1,                'S' => 0 },
            'TTC' => { 'N' => 1.5,              'S' => 0.5 },
            'TTG' => { 'N' => 1.5,              'S' => 0.5 },
            'TTT' => { 'N' => 1.5,              'S' => 0.5 },
        },
        'ATC' => {
            'AAA' => { 'N' => 1.5,              'S' => 0.5 },
            'AAC' => { 'N' => 1,                'S' => 0 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 1 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 0 },
            'ACG' => { 'N' => 1.5,              'S' => 0.5 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 1,                'S' => 0 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 1,                'S' => 1 },
            'ATA' => { 'N' => 0,                'S' => 1 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 0,                'S' => 1 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 0 },
            'CTG' => { 'N' => 1.5,              'S' => 0.5 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2.5,              'S' => 0.5 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 0 },
            'GTG' => { 'N' => 1.5,              'S' => 0.5 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2.5,              'S' => 0.5 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 3,                'S' => 0 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 1.5,              'S' => 0.5 },
            'TTC' => { 'N' => 1,                'S' => 0 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 1,                'S' => 1 },
        },
        'ATG' => {
            'AAA' => { 'N' => 1.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 1,                'S' => 0 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 1.5,              'S' => 0.5 },
            'ACC' => { 'N' => 1.5,              'S' => 0.5 },
            'ACG' => { 'N' => 1,                'S' => 0 },
            'ACT' => { 'N' => 1.5,              'S' => 0.5 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 1,                'S' => 0 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 1,                'S' => 0 },
            'ATC' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 1,                'S' => 0 },
            'CAA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CAC' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'CCA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CCC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CCG' => { 'N' => 2,                'S' => 0 },
            'CCT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGA' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CGC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGG' => { 'N' => 1.5,              'S' => 0.5 },
            'CGT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CTA' => { 'N' => 1.5,              'S' => 0.5 },
            'CTC' => { 'N' => 1.5,              'S' => 0.5 },
            'CTG' => { 'N' => 1,                'S' => 0 },
            'CTT' => { 'N' => 1.5,              'S' => 0.5 },
            'GAA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GAC' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GCA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCG' => { 'N' => 2,                'S' => 0 },
            'GCT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGC' => { 'N' => 2.5,              'S' => 0.5 },
            'GGG' => { 'N' => 2,                'S' => 0 },
            'GGT' => { 'N' => 2.5,              'S' => 0.5 },
            'GTA' => { 'N' => 1.5,              'S' => 0.5 },
            'GTC' => { 'N' => 1.5,              'S' => 0.5 },
            'GTG' => { 'N' => 1,                'S' => 0 },
            'GTT' => { 'N' => 1.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 3,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 3,                'S' => 0 },
            'TCA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 0 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 3,                'S' => 0 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 3,                'S' => 0 },
            'TTA' => { 'N' => 1.5,              'S' => 0.5 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 1,                'S' => 0 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'ATT' => {
            'AAA' => { 'N' => 1.5,              'S' => 0.5 },
            'AAC' => { 'N' => 1,                'S' => 1 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1.5,              'S' => 0.5 },
            'ACT' => { 'N' => 1,                'S' => 0 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 1,                'S' => 1 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 0,                'S' => 1 },
            'ATC' => { 'N' => 0,                'S' => 1 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1.5,              'S' => 0.5 },
            'CTT' => { 'N' => 1,                'S' => 0 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2.5,              'S' => 0.5 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1.5,              'S' => 0.5 },
            'GTT' => { 'N' => 1,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2.5,              'S' => 0.5 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 3,                'S' => 0 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 1.5,              'S' => 0.5 },
            'TTC' => { 'N' => 1,                'S' => 1 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 1,                'S' => 0 },
        },
        'CAA' => {
            'AAA' => { 'N' => 1,                'S' => 0 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 1,                'S' => 1 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2,                'S' => 0 },
            'ACC' => { 'N' => 2.5,              'S' => 0.5 },
            'ACG' => { 'N' => 2,                'S' => 1 },
            'ACT' => { 'N' => 2.5,              'S' => 0.5 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AGG' => { 'N' => 1.5,              'S' => 1.5 },
            'AGT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'ATA' => { 'N' => 2,                'S' => 0 },
            'ATC' => { 'N' => 2.5,              'S' => 0.5 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 1,                'S' => 0 },
            'CAG' => { 'N' => 0,                'S' => 1 },
            'CAT' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 1,                'S' => 0 },
            'CCC' => { 'N' => 1.5,              'S' => 0.5 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1.5,              'S' => 0.5 },
            'CGA' => { 'N' => 1,                'S' => 0 },
            'CGC' => { 'N' => 1.5,              'S' => 0.5 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1.5,              'S' => 0.5 },
            'CTA' => { 'N' => 1,                'S' => 0 },
            'CTC' => { 'N' => 1.5,              'S' => 0.5 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1.5,              'S' => 0.5 },
            'GAA' => { 'N' => 1,                'S' => 0 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 1,                'S' => 1 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2,                'S' => 0 },
            'GCC' => { 'N' => 2.5,              'S' => 0.5 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2.5,              'S' => 0.5 },
            'GGA' => { 'N' => 2,                'S' => 0 },
            'GGC' => { 'N' => 2.5,              'S' => 0.5 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2.5,              'S' => 0.5 },
            'GTA' => { 'N' => 2,                'S' => 0 },
            'GTC' => { 'N' => 2.5,              'S' => 0.5 },
            'GTG' => { 'N' => 2,                'S' => 1 },
            'GTT' => { 'N' => 2.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 0 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'TGG' => { 'N' => 2,                'S' => 1 },
            'TGT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'TTA' => { 'N' => 1,                'S' => 1 },
            'TTC' => { 'N' => 2.5,              'S' => 0.5 },
            'TTG' => { 'N' => 1,                'S' => 2 },
            'TTT' => { 'N' => 2.5,              'S' => 0.5 },
        },
        'CAC' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 1,                'S' => 0 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 1 },
            'ACA' => { 'N' => 2.5,              'S' => 0.5 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2.5,              'S' => 0.5 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 2.5,              'S' => 0.5 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 2.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 1,                'S' => 0 },
            'CAG' => { 'N' => 1,                'S' => 0 },
            'CAT' => { 'N' => 0,                'S' => 1 },
            'CCA' => { 'N' => 1.5,              'S' => 0.5 },
            'CCC' => { 'N' => 1,                'S' => 0 },
            'CCG' => { 'N' => 1.5,              'S' => 0.5 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 1.5,              'S' => 0.5 },
            'CGC' => { 'N' => 1,                'S' => 0 },
            'CGG' => { 'N' => 1.5,              'S' => 0.5 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 1.5,              'S' => 0.5 },
            'CTC' => { 'N' => 1,                'S' => 0 },
            'CTG' => { 'N' => 1.5,              'S' => 0.5 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 1,                'S' => 0 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 1,                'S' => 1 },
            'GCA' => { 'N' => 2.5,              'S' => 0.5 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2.5,              'S' => 0.5 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2.5,              'S' => 0.5 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2.5,              'S' => 0.5 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2.5,              'S' => 0.5 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2.5,              'S' => 0.5 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 1,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 1,                'S' => 1 },
            'TCA' => { 'N' => 2.25,             'S' => 0.75 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2.25,             'S' => 0.75 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 2.75,             'S' => 0.25 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 2.25,             'S' => 0.75 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 2.25,             'S' => 0.75 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'CAG' => {
            'AAA' => { 'N' => 1,                'S' => 1 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 1,                'S' => 0 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2.5,              'S' => 0.5 },
            'ACG' => { 'N' => 2,                'S' => 0 },
            'ACT' => { 'N' => 2.5,              'S' => 0.5 },
            'AGA' => { 'N' => 1.5,              'S' => 1.5 },
            'AGC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AGG' => { 'N' => 1.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'ATA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATC' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'CAA' => { 'N' => 0,                'S' => 1 },
            'CAC' => { 'N' => 1,                'S' => 0 },
            'CAT' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1.5,              'S' => 0.5 },
            'CCG' => { 'N' => 1,                'S' => 0 },
            'CCT' => { 'N' => 1.5,              'S' => 0.5 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1.5,              'S' => 0.5 },
            'CGG' => { 'N' => 1,                'S' => 0 },
            'CGT' => { 'N' => 1.5,              'S' => 0.5 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1.5,              'S' => 0.5 },
            'CTG' => { 'N' => 1,                'S' => 0 },
            'CTT' => { 'N' => 1.5,              'S' => 0.5 },
            'GAA' => { 'N' => 1,                'S' => 1 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 1,                'S' => 0 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2.5,              'S' => 0.5 },
            'GCG' => { 'N' => 2,                'S' => 0 },
            'GCT' => { 'N' => 2.5,              'S' => 0.5 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2.5,              'S' => 0.5 },
            'GGG' => { 'N' => 2,                'S' => 0 },
            'GGT' => { 'N' => 2.5,              'S' => 0.5 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2.5,              'S' => 0.5 },
            'GTG' => { 'N' => 2,                'S' => 0 },
            'GTT' => { 'N' => 2.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 0 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2.75,             'S' => 0.25 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 2.75,             'S' => 0.25 },
            'TTA' => { 'N' => 1,                'S' => 2 },
            'TTC' => { 'N' => 2.5,              'S' => 0.5 },
            'TTG' => { 'N' => 1,                'S' => 1 },
            'TTT' => { 'N' => 2.5,              'S' => 0.5 },
        },
        'CAT' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 1,                'S' => 1 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 2.5,              'S' => 0.5 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.5,              'S' => 0.5 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 2.5,              'S' => 0.5 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 2.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 1,                'S' => 0 },
            'CAC' => { 'N' => 0,                'S' => 1 },
            'CAG' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 1.5,              'S' => 0.5 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1.5,              'S' => 0.5 },
            'CCT' => { 'N' => 1,                'S' => 0 },
            'CGA' => { 'N' => 1.5,              'S' => 0.5 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1.5,              'S' => 0.5 },
            'CGT' => { 'N' => 1,                'S' => 0 },
            'CTA' => { 'N' => 1.5,              'S' => 0.5 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1.5,              'S' => 0.5 },
            'CTT' => { 'N' => 1,                'S' => 0 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 1,                'S' => 1 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 2.5,              'S' => 0.5 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2.5,              'S' => 0.5 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 2.5,              'S' => 0.5 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2.5,              'S' => 0.5 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 2.5,              'S' => 0.5 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2.5,              'S' => 0.5 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 1,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 1,                'S' => 0 },
            'TCA' => { 'N' => 2.25,             'S' => 0.75 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2.25,             'S' => 0.75 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2.75,             'S' => 0.25 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 2.25,             'S' => 0.75 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 2.25,             'S' => 0.75 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'CCA' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 2.5,              'S' => 0.5 },
            'AAG' => { 'N' => 2,                'S' => 1 },
            'AAT' => { 'N' => 2.5,              'S' => 0.5 },
            'ACA' => { 'N' => 1,                'S' => 0 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGG' => { 'N' => 1.5,              'S' => 1.5 },
            'AGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATA' => { 'N' => 2,                'S' => 0 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 1,                'S' => 0 },
            'CAC' => { 'N' => 1.5,              'S' => 0.5 },
            'CAG' => { 'N' => 1,                'S' => 1 },
            'CAT' => { 'N' => 1.5,              'S' => 0.5 },
            'CCC' => { 'N' => 0,                'S' => 1 },
            'CCG' => { 'N' => 0,                'S' => 1 },
            'CCT' => { 'N' => 0,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 0 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 0 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2,                'S' => 1 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 1,                'S' => 0 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 0 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2,                'S' => 0 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2,                'S' => 1 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.25,             'S' => 0.75 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.25,             'S' => 0.75 },
            'TCA' => { 'N' => 1,                'S' => 0 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2,                'S' => 1 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 1.5,              'S' => 0.5 },
            'TTC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TTG' => { 'N' => 1.5,              'S' => 1.5 },
            'TTT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
        },
        'CCC' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 2.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 0 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 1 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 2,                'S' => 1 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 2,                'S' => 1 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 1.5,              'S' => 0.5 },
            'CAC' => { 'N' => 1,                'S' => 0 },
            'CAG' => { 'N' => 1.5,              'S' => 0.5 },
            'CAT' => { 'N' => 1,                'S' => 1 },
            'CCA' => { 'N' => 0,                'S' => 1 },
            'CCG' => { 'N' => 0,                'S' => 1 },
            'CCT' => { 'N' => 0,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 0 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 0 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 0 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2,                'S' => 1 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 0 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 2,                'S' => 1 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 2,                'S' => 1 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'CCG' => {
            'AAA' => { 'N' => 2,                'S' => 1 },
            'AAC' => { 'N' => 2.5,              'S' => 0.5 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 2.5,              'S' => 0.5 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 0 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 1.5,              'S' => 1.5 },
            'AGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGG' => { 'N' => 1.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CAA' => { 'N' => 1,                'S' => 1 },
            'CAC' => { 'N' => 1.5,              'S' => 0.5 },
            'CAG' => { 'N' => 1,                'S' => 0 },
            'CAT' => { 'N' => 1.5,              'S' => 0.5 },
            'CCA' => { 'N' => 0,                'S' => 1 },
            'CCC' => { 'N' => 0,                'S' => 1 },
            'CCT' => { 'N' => 0,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 0 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 0 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 1 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 0 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2,                'S' => 0 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2,                'S' => 0 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.25,             'S' => 0.75 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.25,             'S' => 0.75 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 0 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTA' => { 'N' => 1.5,              'S' => 1.5 },
            'TTC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TTG' => { 'N' => 1.5,              'S' => 0.5 },
            'TTT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
        },
        'CCT' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 2.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 0 },
            'AGA' => { 'N' => 2,                'S' => 1 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 2,                'S' => 1 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 2,                'S' => 1 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 1.5,              'S' => 0.5 },
            'CAC' => { 'N' => 1,                'S' => 1 },
            'CAG' => { 'N' => 1.5,              'S' => 0.5 },
            'CAT' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 0,                'S' => 1 },
            'CCC' => { 'N' => 0,                'S' => 1 },
            'CCG' => { 'N' => 0,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 0 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 0 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 0 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2,                'S' => 1 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 2,                'S' => 1 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 2,                'S' => 1 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'CGA' => {
            'AAA' => { 'N' => 1.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AAG' => { 'N' => 1.5,              'S' => 1.5 },
            'AAT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ACA' => { 'N' => 1.5,              'S' => 0.5 },
            'ACC' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'ACG' => { 'N' => 1.5,              'S' => 1.5 },
            'ACT' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'AGA' => { 'N' => 0,                'S' => 1 },
            'AGC' => { 'N' => 1,                'S' => 1 },
            'AGG' => { 'N' => 0,                'S' => 2 },
            'AGT' => { 'N' => 1,                'S' => 1 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'ATG' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'ATT' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CAA' => { 'N' => 1,                'S' => 0 },
            'CAC' => { 'N' => 1.5,              'S' => 0.5 },
            'CAG' => { 'N' => 1,                'S' => 1 },
            'CAT' => { 'N' => 1.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1,                'S' => 0 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 0,                'S' => 1 },
            'CGG' => { 'N' => 0,                'S' => 1 },
            'CGT' => { 'N' => 0,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 0 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2,                'S' => 1 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 2,                'S' => 0 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 0 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 2,                'S' => 0 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2,                'S' => 1 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TCA' => { 'N' => 2,                'S' => 0 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 1,                'S' => 1 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 1,                'S' => 1 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 1.25,             'S' => 1.75 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'CGC' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 2.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2,                'S' => 1 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 1,                'S' => 1 },
            'AGC' => { 'N' => 1,                'S' => 0 },
            'AGG' => { 'N' => 1,                'S' => 1 },
            'AGT' => { 'N' => 1,                'S' => 1 },
            'ATA' => { 'N' => 2,                'S' => 1 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 1.5,              'S' => 0.5 },
            'CAC' => { 'N' => 1,                'S' => 0 },
            'CAG' => { 'N' => 1.5,              'S' => 0.5 },
            'CAT' => { 'N' => 1,                'S' => 1 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 0 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 0,                'S' => 1 },
            'CGG' => { 'N' => 0,                'S' => 1 },
            'CGT' => { 'N' => 0,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 0 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 0 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2,                'S' => 1 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 0 },
            'TGG' => { 'N' => 1.5,              'S' => 0.5 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 2,                'S' => 1 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'CGG' => {
            'AAA' => { 'N' => 1.5,              'S' => 1.5 },
            'AAC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AAG' => { 'N' => 1.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ACA' => { 'N' => 1.5,              'S' => 1.5 },
            'ACC' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'ACG' => { 'N' => 1.5,              'S' => 0.5 },
            'ACT' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'AGA' => { 'N' => 0,                'S' => 2 },
            'AGC' => { 'N' => 1,                'S' => 1 },
            'AGG' => { 'N' => 0,                'S' => 1 },
            'AGT' => { 'N' => 1,                'S' => 1 },
            'ATA' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'ATC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CAA' => { 'N' => 1,                'S' => 1 },
            'CAC' => { 'N' => 1.5,              'S' => 0.5 },
            'CAG' => { 'N' => 1,                'S' => 0 },
            'CAT' => { 'N' => 1.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 0 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 0,                'S' => 1 },
            'CGC' => { 'N' => 0,                'S' => 1 },
            'CGT' => { 'N' => 0,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 0 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 1 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2,                'S' => 0 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 0 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2,                'S' => 0 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.5,              'S' => 0.5 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.5,              'S' => 0.5 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCG' => { 'N' => 2,                'S' => 0 },
            'TCT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 1.5,              'S' => 0.5 },
            'TGG' => { 'N' => 1,                'S' => 0 },
            'TGT' => { 'N' => 1.5,              'S' => 0.5 },
            'TTA' => { 'N' => 1.25,             'S' => 1.75 },
            'TTC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTG' => { 'N' => 1.5,              'S' => 0.5 },
            'TTT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
        },
        'CGT' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 2.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2,                'S' => 1 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 1,                'S' => 1 },
            'AGC' => { 'N' => 1,                'S' => 1 },
            'AGG' => { 'N' => 1,                'S' => 1 },
            'AGT' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 2,                'S' => 1 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 1.5,              'S' => 0.5 },
            'CAC' => { 'N' => 1,                'S' => 1 },
            'CAG' => { 'N' => 1.5,              'S' => 0.5 },
            'CAT' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 0 },
            'CGA' => { 'N' => 0,                'S' => 1 },
            'CGC' => { 'N' => 0,                'S' => 1 },
            'CGG' => { 'N' => 0,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 0 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 0 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2,                'S' => 1 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 1.5,              'S' => 0.5 },
            'TGT' => { 'N' => 1,                'S' => 0 },
            'TTA' => { 'N' => 2,                'S' => 1 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'CTA' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 2.5,              'S' => 0.5 },
            'AAG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AAT' => { 'N' => 2.5,              'S' => 0.5 },
            'ACA' => { 'N' => 2,                'S' => 0 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGG' => { 'N' => 1.66666666666667, 'S' => 1.33333333333333 },
            'AGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATA' => { 'N' => 1,                'S' => 0 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 1,                'S' => 0 },
            'CAC' => { 'N' => 1.5,              'S' => 0.5 },
            'CAG' => { 'N' => 1,                'S' => 1 },
            'CAT' => { 'N' => 1.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1,                'S' => 0 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 0 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 0,                'S' => 1 },
            'CTG' => { 'N' => 0,                'S' => 1 },
            'CTT' => { 'N' => 0,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2,                'S' => 1 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 2,                'S' => 0 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 0 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 0 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },               # STOP
            'TAC' => { 'N' => 2.25,             'S' => 0.75 },
            'TAG' => { 'N' => 0,                'S' => 0 },               # STOP
            'TAT' => { 'N' => 2.25,             'S' => 0.75 },
            'TCA' => { 'N' => 1.5,              'S' => 0.5 },
            'TCC' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'TCG' => { 'N' => 1.5,              'S' => 1.5 },
            'TCT' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'TGA' => { 'N' => 0,                'S' => 0 },               # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 1.5,              'S' => 1.5 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 0,                'S' => 1 },
            'TTC' => { 'N' => 1,                'S' => 1 },
            'TTG' => { 'N' => 0,                'S' => 2 },
            'TTT' => { 'N' => 1,                'S' => 1 },
        },
        'CTC' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 1 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 1,                'S' => 1 },
            'ATC' => { 'N' => 1,                'S' => 0 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 1.5,              'S' => 0.5 },
            'CAC' => { 'N' => 1,                'S' => 0 },
            'CAG' => { 'N' => 1.5,              'S' => 0.5 },
            'CAT' => { 'N' => 1,                'S' => 1 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 0 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 0 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 0,                'S' => 1 },
            'CTG' => { 'N' => 0,                'S' => 1 },
            'CTT' => { 'N' => 0,                'S' => 1 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 0 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 1,                'S' => 1 },
            'TTC' => { 'N' => 1,                'S' => 0 },
            'TTG' => { 'N' => 1,                'S' => 1 },
            'TTT' => { 'N' => 1,                'S' => 1 },
        },
        'CTG' => {
            'AAA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AAC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'ACA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACG' => { 'N' => 2,                'S' => 0 },
            'ACT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGA' => { 'N' => 1.66666666666667, 'S' => 1.33333333333333 },
            'AGC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGG' => { 'N' => 1.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1.5,              'S' => 0.5 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 1.5,              'S' => 0.5 },
            'CAA' => { 'N' => 1,                'S' => 1 },
            'CAC' => { 'N' => 1.5,              'S' => 0.5 },
            'CAG' => { 'N' => 1,                'S' => 0 },
            'CAT' => { 'N' => 1.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 0 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 0 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 0,                'S' => 1 },
            'CTC' => { 'N' => 0,                'S' => 1 },
            'CTT' => { 'N' => 0,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 1 },
            'GAC' => { 'N' => 2.5,              'S' => 0.5 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.5,              'S' => 0.5 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2,                'S' => 0 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2,                'S' => 0 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 0 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },               # STOP
            'TAC' => { 'N' => 2.25,             'S' => 0.75 },
            'TAG' => { 'N' => 0,                'S' => 0 },               # STOP
            'TAT' => { 'N' => 2.25,             'S' => 0.75 },
            'TCA' => { 'N' => 1.5,              'S' => 1.5 },
            'TCC' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'TCG' => { 'N' => 1.5,              'S' => 0.5 },
            'TCT' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'TGA' => { 'N' => 0,                'S' => 0 },               # STOP
            'TGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TGG' => { 'N' => 1.5,              'S' => 0.5 },
            'TGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TTA' => { 'N' => 0,                'S' => 2 },
            'TTC' => { 'N' => 1,                'S' => 1 },
            'TTG' => { 'N' => 0,                'S' => 1 },
            'TTT' => { 'N' => 1,                'S' => 1 },
        },
        'CTT' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 2,                'S' => 1 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 1,                'S' => 1 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 0 },
            'CAA' => { 'N' => 1.5,              'S' => 0.5 },
            'CAC' => { 'N' => 1,                'S' => 1 },
            'CAG' => { 'N' => 1.5,              'S' => 0.5 },
            'CAT' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 0 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 0 },
            'CTA' => { 'N' => 0,                'S' => 1 },
            'CTC' => { 'N' => 0,                'S' => 1 },
            'CTG' => { 'N' => 0,                'S' => 1 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 1,                'S' => 1 },
            'TTC' => { 'N' => 1,                'S' => 1 },
            'TTG' => { 'N' => 1,                'S' => 1 },
            'TTT' => { 'N' => 1,                'S' => 0 },
        },
        'GAA' => {
            'AAA' => { 'N' => 1,                'S' => 0 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 1,                'S' => 1 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2,                'S' => 0 },
            'ACC' => { 'N' => 2.5,              'S' => 0.5 },
            'ACG' => { 'N' => 2,                'S' => 1 },
            'ACT' => { 'N' => 2.5,              'S' => 0.5 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'AGG' => { 'N' => 2,                'S' => 1 },
            'AGT' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'ATA' => { 'N' => 2,                'S' => 0 },
            'ATC' => { 'N' => 2.5,              'S' => 0.5 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2.5,              'S' => 0.5 },
            'CAA' => { 'N' => 1,                'S' => 0 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 1,                'S' => 1 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2,                'S' => 0 },
            'CCC' => { 'N' => 2.5,              'S' => 0.5 },
            'CCG' => { 'N' => 2,                'S' => 1 },
            'CCT' => { 'N' => 2.5,              'S' => 0.5 },
            'CGA' => { 'N' => 2,                'S' => 0 },
            'CGC' => { 'N' => 2.5,              'S' => 0.5 },
            'CGG' => { 'N' => 2,                'S' => 1 },
            'CGT' => { 'N' => 2.5,              'S' => 0.5 },
            'CTA' => { 'N' => 2,                'S' => 0 },
            'CTC' => { 'N' => 2.5,              'S' => 0.5 },
            'CTG' => { 'N' => 2,                'S' => 1 },
            'CTT' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 1,                'S' => 0 },
            'GAG' => { 'N' => 0,                'S' => 1 },
            'GAT' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 1,                'S' => 0 },
            'GCC' => { 'N' => 1.5,              'S' => 0.5 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1.5,              'S' => 0.5 },
            'GGA' => { 'N' => 1,                'S' => 0 },
            'GGC' => { 'N' => 1.5,              'S' => 0.5 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1.5,              'S' => 0.5 },
            'GTA' => { 'N' => 1,                'S' => 0 },
            'GTC' => { 'N' => 1.5,              'S' => 0.5 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 0 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'TGG' => { 'N' => 2,                'S' => 1 },
            'TGT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 2.75,             'S' => 0.25 },
            'TTG' => { 'N' => 2,                'S' => 1 },
            'TTT' => { 'N' => 2.75,             'S' => 0.25 },
        },
        'GAC' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 1,                'S' => 0 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 1 },
            'ACA' => { 'N' => 2.5,              'S' => 0.5 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2.5,              'S' => 0.5 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 1,                'S' => 0 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 1,                'S' => 1 },
            'CCA' => { 'N' => 2.5,              'S' => 0.5 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2.5,              'S' => 0.5 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 2.5,              'S' => 0.5 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 2.5,              'S' => 0.5 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2.5,              'S' => 0.5 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 2.5,              'S' => 0.5 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 1,                'S' => 0 },
            'GAG' => { 'N' => 1,                'S' => 0 },
            'GAT' => { 'N' => 0,                'S' => 1 },
            'GCA' => { 'N' => 1.5,              'S' => 0.5 },
            'GCC' => { 'N' => 1,                'S' => 0 },
            'GCG' => { 'N' => 1.5,              'S' => 0.5 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 1.5,              'S' => 0.5 },
            'GGC' => { 'N' => 1,                'S' => 0 },
            'GGG' => { 'N' => 1.5,              'S' => 0.5 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 1.5,              'S' => 0.5 },
            'GTC' => { 'N' => 1,                'S' => 0 },
            'GTG' => { 'N' => 1.5,              'S' => 0.5 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 1,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 1,                'S' => 1 },
            'TCA' => { 'N' => 2.25,             'S' => 0.75 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2.25,             'S' => 0.75 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 2.75,             'S' => 0.25 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 2.75,             'S' => 0.25 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 2.75,             'S' => 0.25 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'GAG' => {
            'AAA' => { 'N' => 1,                'S' => 1 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 1,                'S' => 0 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2.5,              'S' => 0.5 },
            'ACG' => { 'N' => 2,                'S' => 0 },
            'ACT' => { 'N' => 2.5,              'S' => 0.5 },
            'AGA' => { 'N' => 2,                'S' => 1 },
            'AGC' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'ATA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATC' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'CAA' => { 'N' => 1,                'S' => 1 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 1,                'S' => 0 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2.5,              'S' => 0.5 },
            'CCG' => { 'N' => 2,                'S' => 0 },
            'CCT' => { 'N' => 2.5,              'S' => 0.5 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2.5,              'S' => 0.5 },
            'CGG' => { 'N' => 2,                'S' => 0 },
            'CGT' => { 'N' => 2.5,              'S' => 0.5 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2.5,              'S' => 0.5 },
            'CTG' => { 'N' => 2,                'S' => 0 },
            'CTT' => { 'N' => 2.5,              'S' => 0.5 },
            'GAA' => { 'N' => 0,                'S' => 1 },
            'GAC' => { 'N' => 1,                'S' => 0 },
            'GAT' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1.5,              'S' => 0.5 },
            'GCG' => { 'N' => 1,                'S' => 0 },
            'GCT' => { 'N' => 1.5,              'S' => 0.5 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1.5,              'S' => 0.5 },
            'GGG' => { 'N' => 1,                'S' => 0 },
            'GGT' => { 'N' => 1.5,              'S' => 0.5 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1.5,              'S' => 0.5 },
            'GTG' => { 'N' => 1,                'S' => 0 },
            'GTT' => { 'N' => 1.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2.5,              'S' => 0.5 },
            'TCG' => { 'N' => 2,                'S' => 0 },
            'TCT' => { 'N' => 2.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2.75,             'S' => 0.25 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 2.75,             'S' => 0.25 },
            'TTA' => { 'N' => 2,                'S' => 1 },
            'TTC' => { 'N' => 2.75,             'S' => 0.25 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 2.75,             'S' => 0.25 },
        },
        'GAT' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 1,                'S' => 1 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 2.5,              'S' => 0.5 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.5,              'S' => 0.5 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.83333333333333, 'S' => 0.166666666666667 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 1,                'S' => 1 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 2.5,              'S' => 0.5 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2.5,              'S' => 0.5 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 2.5,              'S' => 0.5 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2.5,              'S' => 0.5 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 2.5,              'S' => 0.5 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2.5,              'S' => 0.5 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 1,                'S' => 0 },
            'GAC' => { 'N' => 0,                'S' => 1 },
            'GAG' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 1.5,              'S' => 0.5 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1.5,              'S' => 0.5 },
            'GCT' => { 'N' => 1,                'S' => 0 },
            'GGA' => { 'N' => 1.5,              'S' => 0.5 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1.5,              'S' => 0.5 },
            'GGT' => { 'N' => 1,                'S' => 0 },
            'GTA' => { 'N' => 1.5,              'S' => 0.5 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1.5,              'S' => 0.5 },
            'GTT' => { 'N' => 1,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 1,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 1,                'S' => 0 },
            'TCA' => { 'N' => 2.25,             'S' => 0.75 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2.25,             'S' => 0.75 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2.75,             'S' => 0.25 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 2.75,             'S' => 0.25 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 2.75,             'S' => 0.25 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'GCA' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 2.5,              'S' => 0.5 },
            'AAG' => { 'N' => 2,                'S' => 1 },
            'AAT' => { 'N' => 2.5,              'S' => 0.5 },
            'ACA' => { 'N' => 1,                'S' => 0 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGG' => { 'N' => 2,                'S' => 1 },
            'AGT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATA' => { 'N' => 2,                'S' => 0 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2,                'S' => 1 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1,                'S' => 0 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 0 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2,                'S' => 1 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 0 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2,                'S' => 1 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 1,                'S' => 0 },
            'GAC' => { 'N' => 1.5,              'S' => 0.5 },
            'GAG' => { 'N' => 1,                'S' => 1 },
            'GAT' => { 'N' => 1.5,              'S' => 0.5 },
            'GCC' => { 'N' => 0,                'S' => 1 },
            'GCG' => { 'N' => 0,                'S' => 1 },
            'GCT' => { 'N' => 0,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 0 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 0 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.25,             'S' => 0.75 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.25,             'S' => 0.75 },
            'TCA' => { 'N' => 1,                'S' => 0 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2,                'S' => 1 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTG' => { 'N' => 2,                'S' => 1 },
            'TTT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
        },
        'GCC' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 2.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 0 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 2,                'S' => 1 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 0 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 2,                'S' => 1 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 2,                'S' => 1 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 1.5,              'S' => 0.5 },
            'GAC' => { 'N' => 1,                'S' => 0 },
            'GAG' => { 'N' => 1.5,              'S' => 0.5 },
            'GAT' => { 'N' => 1,                'S' => 1 },
            'GCA' => { 'N' => 0,                'S' => 1 },
            'GCG' => { 'N' => 0,                'S' => 1 },
            'GCT' => { 'N' => 0,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 0 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 0 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 0 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'GCG' => {
            'AAA' => { 'N' => 2,                'S' => 1 },
            'AAC' => { 'N' => 2.5,              'S' => 0.5 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 2.5,              'S' => 0.5 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 0 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 1 },
            'AGC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CAA' => { 'N' => 2,                'S' => 1 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 0 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2,                'S' => 0 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2,                'S' => 0 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 1,                'S' => 1 },
            'GAC' => { 'N' => 1.5,              'S' => 0.5 },
            'GAG' => { 'N' => 1,                'S' => 0 },
            'GAT' => { 'N' => 1.5,              'S' => 0.5 },
            'GCA' => { 'N' => 0,                'S' => 1 },
            'GCC' => { 'N' => 0,                'S' => 1 },
            'GCT' => { 'N' => 0,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 0 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 0 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.25,             'S' => 0.75 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.25,             'S' => 0.75 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 0 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTA' => { 'N' => 2,                'S' => 1 },
            'TTC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
        },
        'GCT' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 2.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 0 },
            'AGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 2,                'S' => 1 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 0 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2,                'S' => 1 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2,                'S' => 1 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 1.5,              'S' => 0.5 },
            'GAC' => { 'N' => 1,                'S' => 1 },
            'GAG' => { 'N' => 1.5,              'S' => 0.5 },
            'GAT' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 0,                'S' => 1 },
            'GCC' => { 'N' => 0,                'S' => 1 },
            'GCG' => { 'N' => 0,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 0 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'GGA' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAG' => { 'N' => 2,                'S' => 1 },
            'AAT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'ACA' => { 'N' => 2,                'S' => 0 },
            'ACC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACG' => { 'N' => 2,                'S' => 1 },
            'ACT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGA' => { 'N' => 1,                'S' => 0 },
            'AGC' => { 'N' => 1.5,              'S' => 0.5 },
            'AGG' => { 'N' => 1,                'S' => 1 },
            'AGT' => { 'N' => 1.5,              'S' => 0.5 },
            'ATA' => { 'N' => 2,                'S' => 0 },
            'ATC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2,                'S' => 1 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 2,                'S' => 0 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2,                'S' => 1 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 0 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 0 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2,                'S' => 1 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 1,                'S' => 0 },
            'GAC' => { 'N' => 1.5,              'S' => 0.5 },
            'GAG' => { 'N' => 1,                'S' => 1 },
            'GAT' => { 'N' => 1.5,              'S' => 0.5 },
            'GCA' => { 'N' => 1,                'S' => 0 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 0,                'S' => 1 },
            'GGG' => { 'N' => 0,                'S' => 1 },
            'GGT' => { 'N' => 0,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 0 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'TCA' => { 'N' => 2,                'S' => 0 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 1,                'S' => 1 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 2.25,             'S' => 0.75 },
            'TTG' => { 'N' => 2,                'S' => 1 },
            'TTT' => { 'N' => 2.25,             'S' => 0.75 },
        },
        'GGC' => {
            'AAA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 1,                'S' => 0 },
            'AGG' => { 'N' => 1.5,              'S' => 0.5 },
            'AGT' => { 'N' => 1,                'S' => 1 },
            'ATA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 2.5,              'S' => 0.5 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2,                'S' => 1 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 0 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 2,                'S' => 1 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 1.5,              'S' => 0.5 },
            'GAC' => { 'N' => 1,                'S' => 0 },
            'GAG' => { 'N' => 1.5,              'S' => 0.5 },
            'GAT' => { 'N' => 1,                'S' => 1 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 0 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 0,                'S' => 1 },
            'GGG' => { 'N' => 0,                'S' => 1 },
            'GGT' => { 'N' => 0,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 0 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 0 },
            'TGG' => { 'N' => 1.5,              'S' => 0.5 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 2.5,              'S' => 0.5 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 2.5,              'S' => 0.5 },
            'TTT' => { 'N' => 2,                'S' => 1 },
        },
        'GGG' => {
            'AAA' => { 'N' => 2,                'S' => 1 },
            'AAC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACG' => { 'N' => 2,                'S' => 0 },
            'ACT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGA' => { 'N' => 1,                'S' => 1 },
            'AGC' => { 'N' => 1.5,              'S' => 0.5 },
            'AGG' => { 'N' => 1,                'S' => 0 },
            'AGT' => { 'N' => 1.5,              'S' => 0.5 },
            'ATA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATC' => { 'N' => 2.5,              'S' => 0.5 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 2.5,              'S' => 0.5 },
            'CAA' => { 'N' => 2,                'S' => 1 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2,                'S' => 0 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 0 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2,                'S' => 0 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 1,                'S' => 1 },
            'GAC' => { 'N' => 1.5,              'S' => 0.5 },
            'GAG' => { 'N' => 1,                'S' => 0 },
            'GAT' => { 'N' => 1.5,              'S' => 0.5 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 0 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 0,                'S' => 1 },
            'GGC' => { 'N' => 0,                'S' => 1 },
            'GGT' => { 'N' => 0,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 0 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2.5,              'S' => 0.5 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2.5,              'S' => 0.5 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCG' => { 'N' => 2,                'S' => 0 },
            'TCT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1.5,              'S' => 0.5 },
            'TGG' => { 'N' => 1,                'S' => 0 },
            'TGT' => { 'N' => 1.5,              'S' => 0.5 },
            'TTA' => { 'N' => 2,                'S' => 1 },
            'TTC' => { 'N' => 2.5,              'S' => 0.5 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 2.5,              'S' => 0.5 },
        },
        'GGT' => {
            'AAA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 1.5,              'S' => 0.5 },
            'AGC' => { 'N' => 1,                'S' => 1 },
            'AGG' => { 'N' => 1.5,              'S' => 0.5 },
            'AGT' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.5,              'S' => 0.5 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2,                'S' => 1 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1,                'S' => 1 },
            'CGT' => { 'N' => 1,                'S' => 0 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2,                'S' => 1 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 1.5,              'S' => 0.5 },
            'GAC' => { 'N' => 1,                'S' => 1 },
            'GAG' => { 'N' => 1.5,              'S' => 0.5 },
            'GAT' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 0 },
            'GGA' => { 'N' => 0,                'S' => 1 },
            'GGC' => { 'N' => 0,                'S' => 1 },
            'GGG' => { 'N' => 0,                'S' => 1 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 1.5,              'S' => 0.5 },
            'TGT' => { 'N' => 1,                'S' => 0 },
            'TTA' => { 'N' => 2.5,              'S' => 0.5 },
            'TTC' => { 'N' => 2,                'S' => 1 },
            'TTG' => { 'N' => 2.5,              'S' => 0.5 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'GTA' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 2.5,              'S' => 0.5 },
            'AAG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AAT' => { 'N' => 2.5,              'S' => 0.5 },
            'ACA' => { 'N' => 2,                'S' => 0 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATA' => { 'N' => 1,                'S' => 0 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2,                'S' => 1 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 2,                'S' => 0 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2,                'S' => 1 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 0 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2,                'S' => 1 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 0 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 1,                'S' => 0 },
            'GAC' => { 'N' => 1.5,              'S' => 0.5 },
            'GAG' => { 'N' => 1,                'S' => 1 },
            'GAT' => { 'N' => 1.5,              'S' => 0.5 },
            'GCA' => { 'N' => 1,                'S' => 0 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 0 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 0,                'S' => 1 },
            'GTG' => { 'N' => 0,                'S' => 1 },
            'GTT' => { 'N' => 0,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAC' => { 'N' => 2.5,              'S' => 0.5 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2.5,              'S' => 0.5 },
            'TCA' => { 'N' => 2,                'S' => 0 },
            'TCC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCG' => { 'N' => 2,                'S' => 1 },
            'TCT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 2.25,             'S' => 0.75 },
            'TGG' => { 'N' => 2,                'S' => 1 },
            'TGT' => { 'N' => 2.25,             'S' => 0.75 },
            'TTA' => { 'N' => 1,                'S' => 0 },
            'TTC' => { 'N' => 1.5,              'S' => 0.5 },
            'TTG' => { 'N' => 1,                'S' => 1 },
            'TTT' => { 'N' => 1.5,              'S' => 0.5 },
        },
        'GTC' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 2.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 1,                'S' => 1 },
            'ATC' => { 'N' => 1,                'S' => 0 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2,                'S' => 1 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 2,                'S' => 1 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 0 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 1.5,              'S' => 0.5 },
            'GAC' => { 'N' => 1,                'S' => 0 },
            'GAG' => { 'N' => 1.5,              'S' => 0.5 },
            'GAT' => { 'N' => 1,                'S' => 1 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 0 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 0 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 0,                'S' => 1 },
            'GTG' => { 'N' => 0,                'S' => 1 },
            'GTT' => { 'N' => 0,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 1 },
            'TCA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCC' => { 'N' => 2,                'S' => 0 },
            'TCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCT' => { 'N' => 2,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 2.5,              'S' => 0.5 },
            'TGT' => { 'N' => 2,                'S' => 1 },
            'TTA' => { 'N' => 1.5,              'S' => 0.5 },
            'TTC' => { 'N' => 1,                'S' => 0 },
            'TTG' => { 'N' => 1.5,              'S' => 0.5 },
            'TTT' => { 'N' => 1,                'S' => 1 },
        },
        'GTG' => {
            'AAA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AAC' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'ACA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACG' => { 'N' => 2,                'S' => 0 },
            'ACT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'AGC' => { 'N' => 2.5,              'S' => 0.5 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 2.5,              'S' => 0.5 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1.5,              'S' => 0.5 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 1.5,              'S' => 0.5 },
            'CAA' => { 'N' => 2,                'S' => 1 },
            'CAC' => { 'N' => 2.5,              'S' => 0.5 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 2.5,              'S' => 0.5 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2,                'S' => 0 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2,                'S' => 0 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 0 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 1,                'S' => 1 },
            'GAC' => { 'N' => 1.5,              'S' => 0.5 },
            'GAG' => { 'N' => 1,                'S' => 0 },
            'GAT' => { 'N' => 1.5,              'S' => 0.5 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 0 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 0 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 0,                'S' => 1 },
            'GTC' => { 'N' => 0,                'S' => 1 },
            'GTT' => { 'N' => 0,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2.5,              'S' => 0.5 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2.5,              'S' => 0.5 },
            'TCA' => { 'N' => 2,                'S' => 1 },
            'TCC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCG' => { 'N' => 2,                'S' => 0 },
            'TCT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2.5,              'S' => 0.5 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 2.5,              'S' => 0.5 },
            'TTA' => { 'N' => 1,                'S' => 1 },
            'TTC' => { 'N' => 1.5,              'S' => 0.5 },
            'TTG' => { 'N' => 1,                'S' => 0 },
            'TTT' => { 'N' => 1.5,              'S' => 0.5 },
        },
        'GTT' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 2.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 1,                'S' => 1 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1,                'S' => 0 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2,                'S' => 1 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2,                'S' => 1 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 0 },
            'GAA' => { 'N' => 1.5,              'S' => 0.5 },
            'GAC' => { 'N' => 1,                'S' => 1 },
            'GAG' => { 'N' => 1.5,              'S' => 0.5 },
            'GAT' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 0 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1,                'S' => 1 },
            'GGT' => { 'N' => 1,                'S' => 0 },
            'GTA' => { 'N' => 0,                'S' => 1 },
            'GTC' => { 'N' => 0,                'S' => 1 },
            'GTG' => { 'N' => 0,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCC' => { 'N' => 2,                'S' => 1 },
            'TCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TCT' => { 'N' => 2,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2,                'S' => 1 },
            'TGG' => { 'N' => 2.5,              'S' => 0.5 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 1.5,              'S' => 0.5 },
            'TTC' => { 'N' => 1,                'S' => 1 },
            'TTG' => { 'N' => 1.5,              'S' => 0.5 },
            'TTT' => { 'N' => 1,                'S' => 0 },
        },
        'TAA' => {    # STOP
            'AAA' => { 'N' => 0,   'S' => 0 },      # STOP
            'AAC' => { 'N' => 0,   'S' => 0 },      # STOP
            'AAG' => { 'N' => 0,   'S' => 0 },      # STOP
            'AAT' => { 'N' => 0,   'S' => 0 },      # STOP
            'ACA' => { 'N' => 0,   'S' => 0 },      # STOP
            'ACC' => { 'N' => 0,   'S' => 0 },      # STOP
            'ACG' => { 'N' => 0,   'S' => 0 },      # STOP
            'ACT' => { 'N' => 0,   'S' => 0 },      # STOP
            'AGA' => { 'N' => 0,   'S' => 0 },      # STOP
            'AGC' => { 'N' => 0,   'S' => 0 },      # STOP
            'AGG' => { 'N' => 0,   'S' => 0 },      # STOP
            'AGT' => { 'N' => 0,   'S' => 0 },      # STOP
            'ATA' => { 'N' => 0,   'S' => 0 },      # STOP
            'ATC' => { 'N' => 0,   'S' => 0 },      # STOP
            'ATG' => { 'N' => 0,   'S' => 0 },      # STOP
            'ATT' => { 'N' => 0,   'S' => 0 },      # STOP
            'CAA' => { 'N' => 0,   'S' => 0 },      # STOP
            'CAC' => { 'N' => 0,   'S' => 0 },      # STOP
            'CAG' => { 'N' => 0,   'S' => 0 },      # STOP
            'CAT' => { 'N' => 0,   'S' => 0 },      # STOP
            'CCA' => { 'N' => 0,   'S' => 0 },      # STOP
            'CCC' => { 'N' => 0,   'S' => 0 },      # STOP
            'CCG' => { 'N' => 0,   'S' => 0 },      # STOP
            'CCT' => { 'N' => 0,   'S' => 0 },      # STOP
            'CGA' => { 'N' => 0,   'S' => 0 },      # STOP
            'CGC' => { 'N' => 0,   'S' => 0 },      # STOP
            'CGG' => { 'N' => 0,   'S' => 0 },      # STOP
            'CGT' => { 'N' => 0,   'S' => 0 },      # STOP
            'CTA' => { 'N' => 0,   'S' => 0 },      # STOP
            'CTC' => { 'N' => 0,   'S' => 0 },      # STOP
            'CTG' => { 'N' => 0,   'S' => 0 },      # STOP
            'CTT' => { 'N' => 0,   'S' => 0 },      # STOP
            'GAA' => { 'N' => 0,   'S' => 0 },      # STOP
            'GAC' => { 'N' => 0,   'S' => 0 },      # STOP
            'GAG' => { 'N' => 0,   'S' => 0 },      # STOP
            'GAT' => { 'N' => 0,   'S' => 0 },      # STOP
            'GCA' => { 'N' => 0,   'S' => 0 },      # STOP
            'GCC' => { 'N' => 0,   'S' => 0 },      # STOP
            'GCG' => { 'N' => 0,   'S' => 0 },      # STOP
            'GCT' => { 'N' => 0,   'S' => 0 },      # STOP
            'GGA' => { 'N' => 0,   'S' => 0 },      # STOP
            'GGC' => { 'N' => 0,   'S' => 0 },      # STOP
            'GGG' => { 'N' => 0,   'S' => 0 },      # STOP
            'GGT' => { 'N' => 0,   'S' => 0 },      # STOP
            'GTA' => { 'N' => 0,   'S' => 0 },      # STOP
            'GTC' => { 'N' => 0,   'S' => 0 },      # STOP
            'GTG' => { 'N' => 0,   'S' => 0 },      # STOP
            'GTT' => { 'N' => 0,   'S' => 0 },      # STOP
            'TAC' => { 'N' => 0,   'S' => 0 },      # STOP
            'TAG' => { 'N' => 0,   'S' => 0 },      # STOP
            'TAT' => { 'N' => 0,   'S' => 0 },      # STOP
            'TCA' => { 'N' => 0,   'S' => 0 },      # STOP
            'TCC' => { 'N' => 0,   'S' => 0 },      # STOP
            'TCG' => { 'N' => 0,   'S' => 0 },      # STOP
            'TCT' => { 'N' => 0,   'S' => 0 },      # STOP
            'TGA' => { 'N' => 0,   'S' => 0 },      # STOP
            'TGC' => { 'N' => 0,   'S' => 0 },      # STOP
            'TGG' => { 'N' => '*', 'S' => '*' },    # STOP
            'TGT' => { 'N' => 0,   'S' => 0 },      # STOP
            'TTA' => { 'N' => 0,   'S' => 0 },      # STOP
            'TTC' => { 'N' => 0,   'S' => 0 },      # STOP
            'TTG' => { 'N' => 0,   'S' => 0 },      # STOP
            'TTT' => { 'N' => 0,   'S' => 0 },      # STOP
        },
        'TAC' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 1,                'S' => 0 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 1 },
            'ACA' => { 'N' => 2.25,             'S' => 0.75 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2.25,             'S' => 0.75 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 3,                'S' => 0 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 3,                'S' => 0 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 3,                'S' => 0 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 1,                'S' => 0 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 1,                'S' => 1 },
            'CCA' => { 'N' => 2.25,             'S' => 0.75 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2.25,             'S' => 0.75 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 2.5,              'S' => 0.5 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 2.25,             'S' => 0.75 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 2.25,             'S' => 0.75 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 1,                'S' => 0 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 1,                'S' => 1 },
            'GCA' => { 'N' => 2.25,             'S' => 0.75 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2.25,             'S' => 0.75 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2.5,              'S' => 0.5 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2.5,              'S' => 0.5 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2.5,              'S' => 0.5 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAG' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAT' => { 'N' => 0,                'S' => 1 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 0 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TGC' => { 'N' => 1,                'S' => 0 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 1,                'S' => 0 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 1,                'S' => 1 },
        },
        'TAG' => {    # STOP
            'AAA' => { 'N' => 0, 'S' => 0 },    # STOP
            'AAC' => { 'N' => 0, 'S' => 0 },    # STOP
            'AAG' => { 'N' => 0, 'S' => 0 },    # STOP
            'AAT' => { 'N' => 0, 'S' => 0 },    # STOP
            'ACA' => { 'N' => 0, 'S' => 0 },    # STOP
            'ACC' => { 'N' => 0, 'S' => 0 },    # STOP
            'ACG' => { 'N' => 0, 'S' => 0 },    # STOP
            'ACT' => { 'N' => 0, 'S' => 0 },    # STOP
            'AGA' => { 'N' => 0, 'S' => 0 },    # STOP
            'AGC' => { 'N' => 0, 'S' => 0 },    # STOP
            'AGG' => { 'N' => 0, 'S' => 0 },    # STOP
            'AGT' => { 'N' => 0, 'S' => 0 },    # STOP
            'ATA' => { 'N' => 0, 'S' => 0 },    # STOP
            'ATC' => { 'N' => 0, 'S' => 0 },    # STOP
            'ATG' => { 'N' => 0, 'S' => 0 },    # STOP
            'ATT' => { 'N' => 0, 'S' => 0 },    # STOP
            'CAA' => { 'N' => 0, 'S' => 0 },    # STOP
            'CAC' => { 'N' => 0, 'S' => 0 },    # STOP
            'CAG' => { 'N' => 0, 'S' => 0 },    # STOP
            'CAT' => { 'N' => 0, 'S' => 0 },    # STOP
            'CCA' => { 'N' => 0, 'S' => 0 },    # STOP
            'CCC' => { 'N' => 0, 'S' => 0 },    # STOP
            'CCG' => { 'N' => 0, 'S' => 0 },    # STOP
            'CCT' => { 'N' => 0, 'S' => 0 },    # STOP
            'CGA' => { 'N' => 0, 'S' => 0 },    # STOP
            'CGC' => { 'N' => 0, 'S' => 0 },    # STOP
            'CGG' => { 'N' => 0, 'S' => 0 },    # STOP
            'CGT' => { 'N' => 0, 'S' => 0 },    # STOP
            'CTA' => { 'N' => 0, 'S' => 0 },    # STOP
            'CTC' => { 'N' => 0, 'S' => 0 },    # STOP
            'CTG' => { 'N' => 0, 'S' => 0 },    # STOP
            'CTT' => { 'N' => 0, 'S' => 0 },    # STOP
            'GAA' => { 'N' => 0, 'S' => 0 },    # STOP
            'GAC' => { 'N' => 0, 'S' => 0 },    # STOP
            'GAG' => { 'N' => 0, 'S' => 0 },    # STOP
            'GAT' => { 'N' => 0, 'S' => 0 },    # STOP
            'GCA' => { 'N' => 0, 'S' => 0 },    # STOP
            'GCC' => { 'N' => 0, 'S' => 0 },    # STOP
            'GCG' => { 'N' => 0, 'S' => 0 },    # STOP
            'GCT' => { 'N' => 0, 'S' => 0 },    # STOP
            'GGA' => { 'N' => 0, 'S' => 0 },    # STOP
            'GGC' => { 'N' => 0, 'S' => 0 },    # STOP
            'GGG' => { 'N' => 0, 'S' => 0 },    # STOP
            'GGT' => { 'N' => 0, 'S' => 0 },    # STOP
            'GTA' => { 'N' => 0, 'S' => 0 },    # STOP
            'GTC' => { 'N' => 0, 'S' => 0 },    # STOP
            'GTG' => { 'N' => 0, 'S' => 0 },    # STOP
            'GTT' => { 'N' => 0, 'S' => 0 },    # STOP
            'TAA' => { 'N' => 0, 'S' => 0 },    # STOP
            'TAC' => { 'N' => 0, 'S' => 0 },    # STOP
            'TAT' => { 'N' => 0, 'S' => 0 },    # STOP
            'TCA' => { 'N' => 0, 'S' => 0 },    # STOP
            'TCC' => { 'N' => 0, 'S' => 0 },    # STOP
            'TCG' => { 'N' => 0, 'S' => 0 },    # STOP
            'TCT' => { 'N' => 0, 'S' => 0 },    # STOP
            'TGA' => { 'N' => 0, 'S' => 0 },    # STOP
            'TGC' => { 'N' => 0, 'S' => 0 },    # STOP
            'TGG' => { 'N' => 0, 'S' => 0 },    # STOP
            'TGT' => { 'N' => 0, 'S' => 0 },    # STOP
            'TTA' => { 'N' => 0, 'S' => 0 },    # STOP
            'TTC' => { 'N' => 0, 'S' => 0 },    # STOP
            'TTG' => { 'N' => 0, 'S' => 0 },    # STOP
            'TTT' => { 'N' => 0, 'S' => 0 },    # STOP
        },
        'TAT' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 1,                'S' => 1 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 1,                'S' => 0 },
            'ACA' => { 'N' => 2.25,             'S' => 0.75 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.25,             'S' => 0.75 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 3,                'S' => 0 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 3,                'S' => 0 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 3,                'S' => 0 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 1,                'S' => 1 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 1,                'S' => 0 },
            'CCA' => { 'N' => 2.25,             'S' => 0.75 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2.25,             'S' => 0.75 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2.5,              'S' => 0.5 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 2.25,             'S' => 0.75 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2.25,             'S' => 0.75 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 1,                'S' => 1 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 1,                'S' => 0 },
            'GCA' => { 'N' => 2.25,             'S' => 0.75 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2.25,             'S' => 0.75 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2.5,              'S' => 0.5 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 2.5,              'S' => 0.5 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2.5,              'S' => 0.5 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TAC' => { 'N' => 0,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },    # STOP
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },    # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 1,                'S' => 0 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 1,                'S' => 1 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 1,                'S' => 0 },
        },
        'TCA' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 2.25,             'S' => 0.75 },
            'AAG' => { 'N' => 2,                'S' => 1 },
            'AAT' => { 'N' => 2.25,             'S' => 0.75 },
            'ACA' => { 'N' => 1,                'S' => 0 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 2.25,             'S' => 0.75 },
            'AGG' => { 'N' => 2,                'S' => 1 },
            'AGT' => { 'N' => 2.25,             'S' => 0.75 },
            'ATA' => { 'N' => 2,                'S' => 0 },
            'ATC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CAA' => { 'N' => 2,                'S' => 0 },
            'CAC' => { 'N' => 2.25,             'S' => 0.75 },
            'CAG' => { 'N' => 2,                'S' => 1 },
            'CAT' => { 'N' => 2.25,             'S' => 0.75 },
            'CCA' => { 'N' => 1,                'S' => 0 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 0 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2,                'S' => 1 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 1.5,              'S' => 0.5 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 1.5,              'S' => 1.5 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 2.25,             'S' => 0.75 },
            'GAG' => { 'N' => 2,                'S' => 1 },
            'GAT' => { 'N' => 2.25,             'S' => 0.75 },
            'GCA' => { 'N' => 1,                'S' => 0 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 0 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2,                'S' => 0 },
            'GTC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTG' => { 'N' => 2,                'S' => 1 },
            'GTT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 1,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 0,                'S' => 1 },
            'TCG' => { 'N' => 0,                'S' => 1 },
            'TCT' => { 'N' => 0,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 1,                'S' => 1 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 1,                'S' => 0 },
            'TTC' => { 'N' => 1.5,              'S' => 0.5 },
            'TTG' => { 'N' => 1,                'S' => 1 },
            'TTT' => { 'N' => 1.5,              'S' => 0.5 },
        },
        'TCC' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 2.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 0 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 2.5,              'S' => 0.5 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 2.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 2.5,              'S' => 0.5 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 0 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 0 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 1,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 1,                'S' => 1 },
            'TCA' => { 'N' => 0,                'S' => 1 },
            'TCG' => { 'N' => 0,                'S' => 1 },
            'TCT' => { 'N' => 0,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 0 },
            'TGG' => { 'N' => 1.5,              'S' => 0.5 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 1.5,              'S' => 0.5 },
            'TTC' => { 'N' => 1,                'S' => 0 },
            'TTG' => { 'N' => 1.5,              'S' => 0.5 },
            'TTT' => { 'N' => 1,                'S' => 1 },
        },
        'TCG' => {
            'AAA' => { 'N' => 2,                'S' => 1 },
            'AAC' => { 'N' => 2.25,             'S' => 0.75 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 2.25,             'S' => 0.75 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 0 },
            'ACT' => { 'N' => 1,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 1 },
            'AGC' => { 'N' => 2.5,              'S' => 0.5 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 2.5,              'S' => 0.5 },
            'ATA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ATC' => { 'N' => 2.5,              'S' => 0.5 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 2.5,              'S' => 0.5 },
            'CAA' => { 'N' => 2,                'S' => 1 },
            'CAC' => { 'N' => 2.25,             'S' => 0.75 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 2.25,             'S' => 0.75 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 0 },
            'CCT' => { 'N' => 1,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CGG' => { 'N' => 2,                'S' => 0 },
            'CGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTA' => { 'N' => 1.5,              'S' => 1.5 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 1.5,              'S' => 0.5 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 1 },
            'GAC' => { 'N' => 2.25,             'S' => 0.75 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.25,             'S' => 0.75 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 0 },
            'GCT' => { 'N' => 1,                'S' => 1 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGG' => { 'N' => 2,                'S' => 0 },
            'GGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTG' => { 'N' => 2,                'S' => 0 },
            'GTT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 1,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 1,                'S' => 1 },
            'TCA' => { 'N' => 0,                'S' => 1 },
            'TCC' => { 'N' => 0,                'S' => 1 },
            'TCT' => { 'N' => 0,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1.5,              'S' => 0.5 },
            'TGG' => { 'N' => 1,                'S' => 0 },
            'TGT' => { 'N' => 1.5,              'S' => 0.5 },
            'TTA' => { 'N' => 1,                'S' => 1 },
            'TTC' => { 'N' => 1.5,              'S' => 0.5 },
            'TTG' => { 'N' => 1,                'S' => 0 },
            'TTT' => { 'N' => 1.5,              'S' => 0.5 },
        },
        'TCT' => {
            'AAA' => { 'N' => 2.5,              'S' => 0.5 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 2.5,              'S' => 0.5 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 1,                'S' => 1 },
            'ACC' => { 'N' => 1,                'S' => 1 },
            'ACG' => { 'N' => 1,                'S' => 1 },
            'ACT' => { 'N' => 1,                'S' => 0 },
            'AGA' => { 'N' => 2.5,              'S' => 0.5 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 2.5,              'S' => 0.5 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 2.5,              'S' => 0.5 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 1,                'S' => 1 },
            'CCC' => { 'N' => 1,                'S' => 1 },
            'CCG' => { 'N' => 1,                'S' => 1 },
            'CCT' => { 'N' => 1,                'S' => 0 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 1.83333333333333, 'S' => 1.16666666666667 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 2.5,              'S' => 0.5 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.5,              'S' => 0.5 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 1,                'S' => 1 },
            'GCC' => { 'N' => 1,                'S' => 1 },
            'GCG' => { 'N' => 1,                'S' => 1 },
            'GCT' => { 'N' => 1,                'S' => 0 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 1,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 1,                'S' => 0 },
            'TCA' => { 'N' => 0,                'S' => 1 },
            'TCC' => { 'N' => 0,                'S' => 1 },
            'TCG' => { 'N' => 0,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 1.5,              'S' => 0.5 },
            'TGT' => { 'N' => 1,                'S' => 0 },
            'TTA' => { 'N' => 1.5,              'S' => 0.5 },
            'TTC' => { 'N' => 1,                'S' => 1 },
            'TTG' => { 'N' => 1.5,              'S' => 0.5 },
            'TTT' => { 'N' => 1,                'S' => 0 },
        },
        'TGA' => {    # STOP
            'AAA' => { 'N' => 0, 'S' => 0 },    # STOP
            'AAC' => { 'N' => 0, 'S' => 0 },    # STOP
            'AAG' => { 'N' => 0, 'S' => 0 },    # STOP
            'AAT' => { 'N' => 0, 'S' => 0 },    # STOP
            'ACA' => { 'N' => 0, 'S' => 0 },    # STOP
            'ACC' => { 'N' => 0, 'S' => 0 },    # STOP
            'ACG' => { 'N' => 0, 'S' => 0 },    # STOP
            'ACT' => { 'N' => 0, 'S' => 0 },    # STOP
            'AGA' => { 'N' => 0, 'S' => 0 },    # STOP
            'AGC' => { 'N' => 0, 'S' => 0 },    # STOP
            'AGG' => { 'N' => 0, 'S' => 0 },    # STOP
            'AGT' => { 'N' => 0, 'S' => 0 },    # STOP
            'ATA' => { 'N' => 0, 'S' => 0 },    # STOP
            'ATC' => { 'N' => 0, 'S' => 0 },    # STOP
            'ATG' => { 'N' => 0, 'S' => 0 },    # STOP
            'ATT' => { 'N' => 0, 'S' => 0 },    # STOP
            'CAA' => { 'N' => 0, 'S' => 0 },    # STOP
            'CAC' => { 'N' => 0, 'S' => 0 },    # STOP
            'CAG' => { 'N' => 0, 'S' => 0 },    # STOP
            'CAT' => { 'N' => 0, 'S' => 0 },    # STOP
            'CCA' => { 'N' => 0, 'S' => 0 },    # STOP
            'CCC' => { 'N' => 0, 'S' => 0 },    # STOP
            'CCG' => { 'N' => 0, 'S' => 0 },    # STOP
            'CCT' => { 'N' => 0, 'S' => 0 },    # STOP
            'CGA' => { 'N' => 0, 'S' => 0 },    # STOP
            'CGC' => { 'N' => 0, 'S' => 0 },    # STOP
            'CGG' => { 'N' => 0, 'S' => 0 },    # STOP
            'CGT' => { 'N' => 0, 'S' => 0 },    # STOP
            'CTA' => { 'N' => 0, 'S' => 0 },    # STOP
            'CTC' => { 'N' => 0, 'S' => 0 },    # STOP
            'CTG' => { 'N' => 0, 'S' => 0 },    # STOP
            'CTT' => { 'N' => 0, 'S' => 0 },    # STOP
            'GAA' => { 'N' => 0, 'S' => 0 },    # STOP
            'GAC' => { 'N' => 0, 'S' => 0 },    # STOP
            'GAG' => { 'N' => 0, 'S' => 0 },    # STOP
            'GAT' => { 'N' => 0, 'S' => 0 },    # STOP
            'GCA' => { 'N' => 0, 'S' => 0 },    # STOP
            'GCC' => { 'N' => 0, 'S' => 0 },    # STOP
            'GCG' => { 'N' => 0, 'S' => 0 },    # STOP
            'GCT' => { 'N' => 0, 'S' => 0 },    # STOP
            'GGA' => { 'N' => 0, 'S' => 0 },    # STOP
            'GGC' => { 'N' => 0, 'S' => 0 },    # STOP
            'GGG' => { 'N' => 0, 'S' => 0 },    # STOP
            'GGT' => { 'N' => 0, 'S' => 0 },    # STOP
            'GTA' => { 'N' => 0, 'S' => 0 },    # STOP
            'GTC' => { 'N' => 0, 'S' => 0 },    # STOP
            'GTG' => { 'N' => 0, 'S' => 0 },    # STOP
            'GTT' => { 'N' => 0, 'S' => 0 },    # STOP
            'TAA' => { 'N' => 0, 'S' => 0 },    # STOP
            'TAC' => { 'N' => 0, 'S' => 0 },    # STOP
            'TAG' => { 'N' => 0, 'S' => 0 },    # STOP
            'TAT' => { 'N' => 0, 'S' => 0 },    # STOP
            'TCA' => { 'N' => 0, 'S' => 0 },    # STOP
            'TCC' => { 'N' => 0, 'S' => 0 },    # STOP
            'TCG' => { 'N' => 0, 'S' => 0 },    # STOP
            'TCT' => { 'N' => 0, 'S' => 0 },    # STOP
            'TGC' => { 'N' => 0, 'S' => 0 },    # STOP
            'TGG' => { 'N' => 0, 'S' => 0 },    # STOP
            'TGT' => { 'N' => 0, 'S' => 0 },    # STOP
            'TTA' => { 'N' => 0, 'S' => 0 },    # STOP
            'TTC' => { 'N' => 0, 'S' => 0 },    # STOP
            'TTG' => { 'N' => 0, 'S' => 0 },    # STOP
            'TTT' => { 'N' => 0, 'S' => 0 },    # STOP
        },
        'TGC' => {
            'AAA' => { 'N' => 3,                'S' => 0 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 3,                'S' => 0 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 2.25,             'S' => 0.75 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2.5,              'S' => 0.5 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 1,                'S' => 0 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 1,                'S' => 1 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 3,                'S' => 0 },
            'ATT' => { 'N' => 2,                'S' => 1 },
            'CAA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.75,             'S' => 0.25 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 0 },
            'CGG' => { 'N' => 1.5,              'S' => 0.5 },
            'CGT' => { 'N' => 1,                'S' => 1 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 0 },
            'CTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTT' => { 'N' => 2,                'S' => 1 },
            'GAA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.75,             'S' => 0.25 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 0 },
            'GGG' => { 'N' => 1.5,              'S' => 0.5 },
            'GGT' => { 'N' => 1,                'S' => 1 },
            'GTA' => { 'N' => 2.25,             'S' => 0.75 },
            'GTC' => { 'N' => 2,                'S' => 0 },
            'GTG' => { 'N' => 2.5,              'S' => 0.5 },
            'GTT' => { 'N' => 2,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 1,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 1,                'S' => 1 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 0 },
            'TCG' => { 'N' => 1.5,              'S' => 0.5 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGG' => { 'N' => 1,                'S' => 0 },
            'TGT' => { 'N' => 0,                'S' => 1 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 1,                'S' => 0 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 1,                'S' => 1 },
        },
        'TGG' => {
            'AAA' => { 'N' => 2,                'S' => 1 },
            'AAC' => { 'N' => 3,                'S' => 0 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 3,                'S' => 0 },
            'ACA' => { 'N' => 2,                'S' => 1 },
            'ACC' => { 'N' => 2.5,              'S' => 0.5 },
            'ACG' => { 'N' => 2,                'S' => 0 },
            'ACT' => { 'N' => 2.5,              'S' => 0.5 },
            'AGA' => { 'N' => 1,                'S' => 1 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 1,                'S' => 0 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 3,                'S' => 0 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 3,                'S' => 0 },
            'CAA' => { 'N' => 2,                'S' => 1 },
            'CAC' => { 'N' => 2.75,             'S' => 0.25 },
            'CAG' => { 'N' => 2,                'S' => 0 },
            'CAT' => { 'N' => 2.75,             'S' => 0.25 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CCG' => { 'N' => 2,                'S' => 0 },
            'CCT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1.5,              'S' => 0.5 },
            'CGG' => { 'N' => 1,                'S' => 0 },
            'CGT' => { 'N' => 1.5,              'S' => 0.5 },
            'CTA' => { 'N' => 1.5,              'S' => 1.5 },
            'CTC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CTG' => { 'N' => 1.5,              'S' => 0.5 },
            'CTT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GAA' => { 'N' => 2,                'S' => 1 },
            'GAC' => { 'N' => 2.75,             'S' => 0.25 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.75,             'S' => 0.25 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCG' => { 'N' => 2,                'S' => 0 },
            'GCT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1.5,              'S' => 0.5 },
            'GGG' => { 'N' => 1,                'S' => 0 },
            'GGT' => { 'N' => 1.5,              'S' => 0.5 },
            'GTA' => { 'N' => 2,                'S' => 1 },
            'GTC' => { 'N' => 2.5,              'S' => 0.5 },
            'GTG' => { 'N' => 2,                'S' => 0 },
            'GTT' => { 'N' => 2.5,              'S' => 0.5 },
            'TAA' => { 'N' => '*',              'S' => '*' },    # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },      # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1.5,              'S' => 0.5 },
            'TCG' => { 'N' => 1,                'S' => 0 },
            'TCT' => { 'N' => 1.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },      # STOP
            'TGC' => { 'N' => 1,                'S' => 0 },
            'TGT' => { 'N' => 1,                'S' => 0 },
            'TTA' => { 'N' => 1,                'S' => 1 },
            'TTC' => { 'N' => 2,                'S' => 0 },
            'TTG' => { 'N' => 1,                'S' => 0 },
            'TTT' => { 'N' => 2,                'S' => 0 },
        },
        'TGT' => {
            'AAA' => { 'N' => 3,                'S' => 0 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 3,                'S' => 0 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2.25,             'S' => 0.75 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.5,              'S' => 0.5 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 1,                'S' => 1 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 1,                'S' => 0 },
            'ATA' => { 'N' => 2.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 1 },
            'ATG' => { 'N' => 3,                'S' => 0 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.75,             'S' => 0.25 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2,                'S' => 1 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 1,                'S' => 1 },
            'CGG' => { 'N' => 1.5,              'S' => 0.5 },
            'CGT' => { 'N' => 1,                'S' => 0 },
            'CTA' => { 'N' => 2,                'S' => 1 },
            'CTC' => { 'N' => 2,                'S' => 1 },
            'CTG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTT' => { 'N' => 2,                'S' => 0 },
            'GAA' => { 'N' => 2.66666666666667, 'S' => 0.333333333333333 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.75,             'S' => 0.25 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 1,                'S' => 1 },
            'GGC' => { 'N' => 1,                'S' => 1 },
            'GGG' => { 'N' => 1.5,              'S' => 0.5 },
            'GGT' => { 'N' => 1,                'S' => 0 },
            'GTA' => { 'N' => 2.25,             'S' => 0.75 },
            'GTC' => { 'N' => 2,                'S' => 1 },
            'GTG' => { 'N' => 2.5,              'S' => 0.5 },
            'GTT' => { 'N' => 2,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 1,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 1,                'S' => 0 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1.5,              'S' => 0.5 },
            'TCT' => { 'N' => 1,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 0,                'S' => 1 },
            'TGG' => { 'N' => 1,                'S' => 0 },
            'TTA' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 1,                'S' => 1 },
            'TTG' => { 'N' => 2,                'S' => 0 },
            'TTT' => { 'N' => 1,                'S' => 0 },
        },
        'TTA' => {
            'AAA' => { 'N' => 2,                'S' => 0 },
            'AAC' => { 'N' => 2.75,             'S' => 0.25 },
            'AAG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AAT' => { 'N' => 2.75,             'S' => 0.25 },
            'ACA' => { 'N' => 2,                'S' => 0 },
            'ACC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ACG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AGA' => { 'N' => 2,                'S' => 0 },
            'AGC' => { 'N' => 2.75,             'S' => 0.25 },
            'AGG' => { 'N' => 2.25,             'S' => 0.75 },
            'AGT' => { 'N' => 2.75,             'S' => 0.25 },
            'ATA' => { 'N' => 1,                'S' => 0 },
            'ATC' => { 'N' => 1.5,              'S' => 0.5 },
            'ATG' => { 'N' => 1.5,              'S' => 0.5 },
            'ATT' => { 'N' => 1.5,              'S' => 0.5 },
            'CAA' => { 'N' => 1,                'S' => 1 },
            'CAC' => { 'N' => 2.25,             'S' => 0.75 },
            'CAG' => { 'N' => 1,                'S' => 2 },
            'CAT' => { 'N' => 2.25,             'S' => 0.75 },
            'CCA' => { 'N' => 1.5,              'S' => 0.5 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 1.5,              'S' => 1.5 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 1.25,             'S' => 1.75 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 0,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 0,                'S' => 2 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 0 },
            'GAC' => { 'N' => 2.75,             'S' => 0.25 },
            'GAG' => { 'N' => 2,                'S' => 1 },
            'GAT' => { 'N' => 2.75,             'S' => 0.25 },
            'GCA' => { 'N' => 2,                'S' => 0 },
            'GCC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCG' => { 'N' => 2,                'S' => 1 },
            'GCT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGA' => { 'N' => 2,                'S' => 0 },
            'GGC' => { 'N' => 2.5,              'S' => 0.5 },
            'GGG' => { 'N' => 2,                'S' => 1 },
            'GGT' => { 'N' => 2.5,              'S' => 0.5 },
            'GTA' => { 'N' => 1,                'S' => 0 },
            'GTC' => { 'N' => 1.5,              'S' => 0.5 },
            'GTG' => { 'N' => 1,                'S' => 1 },
            'GTT' => { 'N' => 1.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 1,                'S' => 0 },
            'TCC' => { 'N' => 1.5,              'S' => 0.5 },
            'TCG' => { 'N' => 1,                'S' => 1 },
            'TCT' => { 'N' => 1.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 1,                'S' => 1 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTC' => { 'N' => 1,                'S' => 0 },
            'TTG' => { 'N' => 0,                'S' => 1 },
            'TTT' => { 'N' => 1,                'S' => 0 },
        },
        'TTC' => {
            'AAA' => { 'N' => 2.75,             'S' => 0.25 },
            'AAC' => { 'N' => 2,                'S' => 0 },
            'AAG' => { 'N' => 3,                'S' => 0 },
            'AAT' => { 'N' => 2,                'S' => 1 },
            'ACA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ACC' => { 'N' => 2,                'S' => 0 },
            'ACG' => { 'N' => 2.5,              'S' => 0.5 },
            'ACT' => { 'N' => 2,                'S' => 1 },
            'AGA' => { 'N' => 2.75,             'S' => 0.25 },
            'AGC' => { 'N' => 2,                'S' => 0 },
            'AGG' => { 'N' => 3,                'S' => 0 },
            'AGT' => { 'N' => 2,                'S' => 1 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1,                'S' => 0 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 1,                'S' => 1 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 0 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 1 },
            'CCA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CCC' => { 'N' => 2,                'S' => 0 },
            'CCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 0 },
            'CGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGT' => { 'N' => 2,                'S' => 1 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 0 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2.75,             'S' => 0.25 },
            'GAC' => { 'N' => 2,                'S' => 0 },
            'GAG' => { 'N' => 2.75,             'S' => 0.25 },
            'GAT' => { 'N' => 2,                'S' => 1 },
            'GCA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCC' => { 'N' => 2,                'S' => 0 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 1 },
            'GGA' => { 'N' => 2.25,             'S' => 0.75 },
            'GGC' => { 'N' => 2,                'S' => 0 },
            'GGG' => { 'N' => 2.5,              'S' => 0.5 },
            'GGT' => { 'N' => 2,                'S' => 1 },
            'GTA' => { 'N' => 1.5,              'S' => 0.5 },
            'GTC' => { 'N' => 1,                'S' => 0 },
            'GTG' => { 'N' => 1.5,              'S' => 0.5 },
            'GTT' => { 'N' => 1,                'S' => 1 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 1,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 1,                'S' => 1 },
            'TCA' => { 'N' => 1.5,              'S' => 0.5 },
            'TCC' => { 'N' => 1,                'S' => 0 },
            'TCG' => { 'N' => 1.5,              'S' => 0.5 },
            'TCT' => { 'N' => 1,                'S' => 1 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 0 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 1,                'S' => 1 },
            'TTA' => { 'N' => 1,                'S' => 0 },
            'TTG' => { 'N' => 1,                'S' => 0 },
            'TTT' => { 'N' => 0,                'S' => 1 },
        },
        'TTG' => {
            'AAA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'AAC' => { 'N' => 3,                'S' => 0 },
            'AAG' => { 'N' => 2,                'S' => 0 },
            'AAT' => { 'N' => 3,                'S' => 0 },
            'ACA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'ACC' => { 'N' => 2.5,              'S' => 0.5 },
            'ACG' => { 'N' => 2,                'S' => 0 },
            'ACT' => { 'N' => 2.5,              'S' => 0.5 },
            'AGA' => { 'N' => 2.25,             'S' => 0.75 },
            'AGC' => { 'N' => 3,                'S' => 0 },
            'AGG' => { 'N' => 2,                'S' => 0 },
            'AGT' => { 'N' => 3,                'S' => 0 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 2,                'S' => 0 },
            'ATG' => { 'N' => 1,                'S' => 0 },
            'ATT' => { 'N' => 2,                'S' => 0 },
            'CAA' => { 'N' => 1,                'S' => 2 },
            'CAC' => { 'N' => 2.25,             'S' => 0.75 },
            'CAG' => { 'N' => 1,                'S' => 1 },
            'CAT' => { 'N' => 2.25,             'S' => 0.75 },
            'CCA' => { 'N' => 1.5,              'S' => 1.5 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 1.5,              'S' => 0.5 },
            'CCT' => { 'N' => 2,                'S' => 1 },
            'CGA' => { 'N' => 1.25,             'S' => 1.75 },
            'CGC' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CGG' => { 'N' => 1.5,              'S' => 0.5 },
            'CGT' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CTA' => { 'N' => 0,                'S' => 2 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 0,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 1 },
            'GAA' => { 'N' => 2,                'S' => 1 },
            'GAC' => { 'N' => 2.75,             'S' => 0.25 },
            'GAG' => { 'N' => 2,                'S' => 0 },
            'GAT' => { 'N' => 2.75,             'S' => 0.25 },
            'GCA' => { 'N' => 2,                'S' => 1 },
            'GCC' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCG' => { 'N' => 2,                'S' => 0 },
            'GCT' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GGA' => { 'N' => 2,                'S' => 1 },
            'GGC' => { 'N' => 2.5,              'S' => 0.5 },
            'GGG' => { 'N' => 2,                'S' => 0 },
            'GGT' => { 'N' => 2.5,              'S' => 0.5 },
            'GTA' => { 'N' => 1,                'S' => 1 },
            'GTC' => { 'N' => 1.5,              'S' => 0.5 },
            'GTG' => { 'N' => 1,                'S' => 0 },
            'GTT' => { 'N' => 1.5,              'S' => 0.5 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 2,                'S' => 0 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 2,                'S' => 0 },
            'TCA' => { 'N' => 1,                'S' => 1 },
            'TCC' => { 'N' => 1.5,              'S' => 0.5 },
            'TCG' => { 'N' => 1,                'S' => 0 },
            'TCT' => { 'N' => 1.5,              'S' => 0.5 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 2,                'S' => 0 },
            'TGG' => { 'N' => 1,                'S' => 0 },
            'TGT' => { 'N' => 2,                'S' => 0 },
            'TTA' => { 'N' => 0,                'S' => 1 },
            'TTC' => { 'N' => 1,                'S' => 0 },
            'TTT' => { 'N' => 1,                'S' => 0 },
        },
        'TTT' => {
            'AAA' => { 'N' => 2.75,             'S' => 0.25 },
            'AAC' => { 'N' => 2,                'S' => 1 },
            'AAG' => { 'N' => 3,                'S' => 0 },
            'AAT' => { 'N' => 2,                'S' => 0 },
            'ACA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'ACC' => { 'N' => 2,                'S' => 1 },
            'ACG' => { 'N' => 2.5,              'S' => 0.5 },
            'ACT' => { 'N' => 2,                'S' => 0 },
            'AGA' => { 'N' => 2.75,             'S' => 0.25 },
            'AGC' => { 'N' => 2,                'S' => 1 },
            'AGG' => { 'N' => 3,                'S' => 0 },
            'AGT' => { 'N' => 2,                'S' => 0 },
            'ATA' => { 'N' => 1.5,              'S' => 0.5 },
            'ATC' => { 'N' => 1,                'S' => 1 },
            'ATG' => { 'N' => 2,                'S' => 0 },
            'ATT' => { 'N' => 1,                'S' => 0 },
            'CAA' => { 'N' => 2.5,              'S' => 0.5 },
            'CAC' => { 'N' => 2,                'S' => 1 },
            'CAG' => { 'N' => 2.5,              'S' => 0.5 },
            'CAT' => { 'N' => 2,                'S' => 0 },
            'CCA' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CCC' => { 'N' => 2,                'S' => 1 },
            'CCG' => { 'N' => 2.16666666666667, 'S' => 0.833333333333333 },
            'CCT' => { 'N' => 2,                'S' => 0 },
            'CGA' => { 'N' => 2,                'S' => 1 },
            'CGC' => { 'N' => 2,                'S' => 1 },
            'CGG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'CGT' => { 'N' => 2,                'S' => 0 },
            'CTA' => { 'N' => 1,                'S' => 1 },
            'CTC' => { 'N' => 1,                'S' => 1 },
            'CTG' => { 'N' => 1,                'S' => 1 },
            'CTT' => { 'N' => 1,                'S' => 0 },
            'GAA' => { 'N' => 2.75,             'S' => 0.25 },
            'GAC' => { 'N' => 2,                'S' => 1 },
            'GAG' => { 'N' => 2.75,             'S' => 0.25 },
            'GAT' => { 'N' => 2,                'S' => 0 },
            'GCA' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCC' => { 'N' => 2,                'S' => 1 },
            'GCG' => { 'N' => 2.33333333333333, 'S' => 0.666666666666667 },
            'GCT' => { 'N' => 2,                'S' => 0 },
            'GGA' => { 'N' => 2.25,             'S' => 0.75 },
            'GGC' => { 'N' => 2,                'S' => 1 },
            'GGG' => { 'N' => 2.5,              'S' => 0.5 },
            'GGT' => { 'N' => 2,                'S' => 0 },
            'GTA' => { 'N' => 1.5,              'S' => 0.5 },
            'GTC' => { 'N' => 1,                'S' => 1 },
            'GTG' => { 'N' => 1.5,              'S' => 0.5 },
            'GTT' => { 'N' => 1,                'S' => 0 },
            'TAA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAC' => { 'N' => 1,                'S' => 1 },
            'TAG' => { 'N' => 0,                'S' => 0 },     # STOP
            'TAT' => { 'N' => 1,                'S' => 0 },
            'TCA' => { 'N' => 1.5,              'S' => 0.5 },
            'TCC' => { 'N' => 1,                'S' => 1 },
            'TCG' => { 'N' => 1.5,              'S' => 0.5 },
            'TCT' => { 'N' => 1,                'S' => 0 },
            'TGA' => { 'N' => 0,                'S' => 0 },     # STOP
            'TGC' => { 'N' => 1,                'S' => 1 },
            'TGG' => { 'N' => 2,                'S' => 0 },
            'TGT' => { 'N' => 1,                'S' => 0 },
            'TTA' => { 'N' => 1,                'S' => 0 },
            'TTC' => { 'N' => 0,                'S' => 1 },
            'TTG' => { 'N' => 1,                'S' => 0 }
        }
    );

    my $num_N_diffs = $all_diffs_hh{$codon1}->{$codon2}->{N};
    my $num_S_diffs = $all_diffs_hh{$codon1}->{$codon2}->{S};

    my @diffs_arr = ( $num_N_diffs, $num_S_diffs );

    return @diffs_arr;
}

#########################################################################################
# End the program by notifying the screen at command line
sub end_the_program {
    my $time2       = time;
    my $local_time2 = localtime;

    my $time_diff              = ( $time2 - $time1 );
    my $time_diff_rounded      = sprintf( "%.2f", $time_diff );
    my $mins_elapsed           = ( $time_diff / 60 );
    my $whole_mins_elapsed     = int($mins_elapsed);
    my $whole_mins_in_secs     = ( $whole_mins_elapsed * 60 );
    my $secs_remaining         = ( $time_diff - $whole_mins_in_secs );
    my $secs_remaining_rounded = sprintf( "%.2f", $secs_remaining );

    print
"SNPGenie completed at local time $local_time2. The process took $time_diff_rounded secs, i.e., "
      . "$whole_mins_elapsed mins and $secs_remaining_rounded secs\n";

    print
"\n################################################################################"
      . "\n##              SNPGenie for Phylogeneous completed successfully.             ##"
      . "\n##        Please find results in the \/data\/ dir and subdirs (families).       ##\n"
      . "################################################################################"
      . "\n\n\n";
}
