#! /share/apps/perl/5.24.0/intel/bin/perl 

#AMNH-HUXLEY-PBS is #! /usr/local/software/PERL/perl-5.26.0/bin/perl
#NYU-PRINCE-SLURM is #! /share/apps/perl/5.24.0/intel/bin/perl
#PC is #! /usr/bin/perl

# CALL FORMAT:
# $OID_HOME/bin/snpgenie_oid_recalc_sw.pl <sw_size> <sw_step_size> <num_bootstraps> <procs_per_node>

# EXAMPLE CALL:
# ./snpgenie_oid_recalc_sw.pl --sliding_window_size=8 --sliding_window_step=1 --num_bootstraps=1000 --procs_per_node=4

# PROGRAM: Recalculate sliding window values using a given window and step size for all families/partitions, for OrthologID and HPC

# Copyright (C) 2017 Chase W. Nelson

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

# DATE CREATED: November 9, 2017

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
# CITATION3: Nelson CW, Hughes AL (2015) Within-host nucleotide diversity of virus
#	populations: Insights from next-generation sequencing. Infection, Genetics and
#	Evolution 30:1-7, doi: 10.1016/j.meegid.2014.11.026

use strict;

#use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw(max);
use Parallel::ForkManager;

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

#########################################################################################
# INITIALIZE VARIABLES (all optional and/or defaulted)
my $sliding_window_size;
my $sliding_window_step;
my $num_bootstraps;
my $procs_per_node;

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions(
    "sliding_window_size=i" => \$sliding_window_size,
    "sliding_window_step=i" => \$sliding_window_step,
    "num_bootstraps=i"      => \$num_bootstraps,
    "procs_per_node:i"      => \$procs_per_node
  )

  or die
"\n### WARNING: Error in command line arguments. SNPGenie for OID terminated.\n\n";

# Get the time
my $time1       = time;
my $local_time1 = localtime;

print "\nSNPGenie sliding windows initiated at local time $local_time1\n\n";

# Test and/or reset OPTIONS given the user's INPUT and RECORD PARAMETERS
if ( !$sliding_window_size ) {
    if ( $sliding_window_size != 0 ) {    # Called as a flag, but given no value
        $sliding_window_size = 10;        # default behavior: 10-mer peptide
    }
}
elsif ( $sliding_window_size < 1 ) {
    die
      "\n### WARNING: The --sliding_window_size option must be an integer ≥1\n"
      . "### SNPGenie for OID terminated.\n\n";
}

if ( !$sliding_window_step ) {
    if ( $sliding_window_step != 0 ) {    # Called as a flag, but given no value
        $sliding_window_step = 1;         # default behavior: one codon
    }
}
elsif ( $sliding_window_step < 1 ) {
    die
      "\n### WARNING: The --sliding_window_step option must be an integer ≥1\n"
      . "### SNPGenie for OID terminated.\n\n";
}

if ( !$num_bootstraps ) {
    $num_bootstraps = 0;    # DEFAULT is 1000
    print "### BOOTSTRAP VALUE NOT GIVEN. DEFAULT IS NO BOOTSTRAPPING.\n"
      . "### (You sure you don't want our statistics, boo?)\n\n";
}
elsif ( $num_bootstraps < 2 ) {
    die "\n### WARNING: The --num_bootstraps option must be an integer ≥2\n"
      . "### SNPGenie for OID terminated.\n\n";
}

if ( !$procs_per_node ) {
    $procs_per_node = 1;    # DEFAULT is 1
    print "### PROCS PER NODE VALUE NOT GIVEN. DEFAULT USED: JUST 1.\n"
      . "### (You sure you don't want our parallelization, boo?)\n\n";
}
elsif ( $procs_per_node < 1 ) {
    die "\n### WARNING: The --procs_per_node option must be an integer ≥1\n"
      . "### SNPGenie for OID terminated.\n\n";
}

#########################################################################################
# GET ALL FAMILY DIRECTORIES (numbers only; no 'S')
chdir("$OID_USER_DIR\/data");
my @OID_USER_DIR_contents = glob "*";
my @OID_USER_DIR_contents_NUMERICAL;
foreach my $family_id (@OID_USER_DIR_contents) {

    #print "\nfamily_id is $family_id\n";
    if ( $family_id =~ /^\d/ )
    {    # only directory names beginning with numbers, please
        push( @OID_USER_DIR_contents_NUMERICAL, $family_id );
    }
}
@OID_USER_DIR_contents = @OID_USER_DIR_contents_NUMERICAL;

#STDOUT->autoflush(1);

#########################################################################################
# Input and output file names; same for every family directory
my $between_group_codon_input_file = 'between_group_codon_results.txt';
my $between_group_sw_output_file =
  'between_group_sw_' . $sliding_window_size . 'codons_results.txt';

#########################################################################################
# IN PARALLEL: LOOP through each family, i.e. directory in /data/
my $pm_family = Parallel::ForkManager->new($procs_per_node);
FAMILY: foreach my $family_id (@OID_USER_DIR_contents) {
    $pm_family->start and next;

    #	my $family_start_time = time;

    print
"#########################################################################################";
    print "\n##### Working on family $family_id\n";
    print
"#########################################################################################";
    print "\n";

    # MAKE SURE WE DON'T OVERWRITE
##	if(-f "$OID_USER_DIR/data/$family_id/$between_group_sw_output_file") { # don't re-do if done
##		print "Family $family_id already completed for this window size. Moving on to next family.\n";
##		#next FAMILY; # this doesn't work because a child process apparently cannot be called
##	} els

    if ( -f "$OID_USER_DIR/data/$family_id/between_group_codon_results.txt" )
    {    # DO IT
         # Go to the family directory
        chdir("$OID_USER_DIR\/data\/$family_id");

        ##################################################################################
        # READY THE OUTPUT FILE

        open( OUTFILE_BETWEEN_SW, ">$between_group_sw_output_file" );
        print OUTFILE_BETWEEN_SW
"analysis\tfamily\tpartition\tgroup_1\tgroup_2\twindow\tfirst_codon\tlast_codon\t"
          . "num_defined_codons_g1\tnum_defined_codons_g2\tmean_num_defined_codons\tmin_num_defined_codons\t"
          . "N_sites\tS_sites\t"
          . "N_diffs\tS_diffs\t"
          . "dN\tdS\t"
          . "dN\-dS\t"
          . "dN\/dS\tdN\>total_dS";

        if ( $num_bootstraps > 1 ) {
            print OUTFILE_BETWEEN_SW "\tSE(dN-dS)\tZ_value\tsignificance\n";
        }
        else {
            print OUTFILE_BETWEEN_SW "\n";
        }

        close OUTFILE_BETWEEN_SW;

        ##################################################################################
        # STORE BETWEEN-GROUP CODON DATA

        # Get header names
        my @header_names_arr =
          &get_header_names( $between_group_codon_input_file,
            $between_group_codon_input_file );

        #print "@header_names_arr";

        my %header_indices;

        # Determine the index of each column
        for ( my $i = 0 ; $i < scalar(@header_names_arr) ; $i++ ) {
            my $curr_header = $header_names_arr[$i];
            $header_indices{$curr_header} = $i;
        }

        # Make sure we got all the expected headers
        my @required_between_group_headers =
          qw/analysis family partition group_1 group_2 codon variability num_defined_codons_g1 num_defined_codons_g2 comparisons N_sites S_sites N_diffs S_diffs/;

        foreach (@required_between_group_headers) {
            unless ( exists $header_indices{$_} ) {
                die
"### DIE: the header name $_ is not present in the between-group codon results file.\n\n";
            }
        }

        # STORE ALL FAMILY BETWEEN-GROUP CODON DATA
        my %between_group_partition_codon_data_hh;
        my %between_group_a_f_p_groups;

        # STORE overall partition values
        my %between_group_a_f_p_sums;

        #		my $between_product_N_sites_sum = 0;
        my $between_product_S_sites_sum = 0;

        #		my $between_product_N_diffs_sum = 0;
        my $between_product_S_diffs_sum = 0;

        open( IN_BETWEEN_GROUP_CODON_FILE, "$between_group_codon_input_file" )
          or die
          "\nCould not open the file $between_group_codon_input_file - $!\n\n";
        my $line = 0;
        while (<IN_BETWEEN_GROUP_CODON_FILE>) {    # record
            if ( $line == 0 ) {
                $line++;
            }
            else {
                chomp;
                my @line_arr = split( /\t/, $_, -1 );

                # Store data here
                my $analysis  = $line_arr[ $header_indices{'analysis'} ];
                my $family    = $line_arr[ $header_indices{'family'} ];
                my $partition = $line_arr[ $header_indices{'partition'} ];
                my $codon     = $line_arr[ $header_indices{'codon'} ];
                $codon =~ s/codon_//;    # in case it's the old format
                $codon =~ s/\s//g;       # in case it's the old format

                $between_group_partition_codon_data_hh{$analysis}->{$family}
                  ->{$partition}->{$codon}->{variability} =
                  $line_arr[ $header_indices{'variability'} ];
                $between_group_partition_codon_data_hh{$analysis}->{$family}
                  ->{$partition}->{$codon}->{num_defined_codons_g1} =
                  $line_arr[ $header_indices{'num_defined_codons_g1'} ];
                $between_group_partition_codon_data_hh{$analysis}->{$family}
                  ->{$partition}->{$codon}->{num_defined_codons_g2} =
                  $line_arr[ $header_indices{'num_defined_codons_g2'} ];
                $between_group_partition_codon_data_hh{$analysis}->{$family}
                  ->{$partition}->{$codon}->{comparisons} =
                  $line_arr[ $header_indices{'comparisons'} ];
                $between_group_partition_codon_data_hh{$analysis}->{$family}
                  ->{$partition}->{$codon}->{N_sites} =
                  $line_arr[ $header_indices{'N_sites'} ];
                $between_group_partition_codon_data_hh{$analysis}->{$family}
                  ->{$partition}->{$codon}->{S_sites} =
                  $line_arr[ $header_indices{'S_sites'} ];
                $between_group_partition_codon_data_hh{$analysis}->{$family}
                  ->{$partition}->{$codon}->{N_diffs} =
                  $line_arr[ $header_indices{'N_diffs'} ];
                $between_group_partition_codon_data_hh{$analysis}->{$family}
                  ->{$partition}->{$codon}->{S_diffs} =
                  $line_arr[ $header_indices{'S_diffs'} ];

#				print "\n\nVariability: " . $between_group_partition_codon_data_hh{$analysis}->{$family}->{$partition}->{$codon}->{variability} . "\n\n";
#				print "\n\nN_diffs: " . $between_group_partition_codon_data_hh{$analysis}->{$family}->{$partition}->{$codon}->{N_diffs} . "\n\n";

#				print "\nSTORING:\nanalysis: $analysis\nfamily: $family\npartition: $partition\ncodon: $codon\nNsites: ".
#					$line_arr[$header_indices{'N_sites'}] . "\n\n";

                # Store group names here
                $between_group_a_f_p_groups{$analysis}->{$family}->{$partition}
                  ->{group_1} = $line_arr[ $header_indices{'group_1'} ];
                $between_group_a_f_p_groups{$analysis}->{$family}->{$partition}
                  ->{group_2} = $line_arr[ $header_indices{'group_2'} ];

                # Sum partition values here
                $between_group_a_f_p_sums{$analysis}->{$family}->{$partition}
                  ->{N_sites} += $line_arr[ $header_indices{'N_sites'} ];
                $between_group_a_f_p_sums{$analysis}->{$family}->{$partition}
                  ->{S_sites} += $line_arr[ $header_indices{'S_sites'} ];
                $between_group_a_f_p_sums{$analysis}->{$family}->{$partition}
                  ->{N_diffs} += $line_arr[ $header_indices{'N_diffs'} ];
                $between_group_a_f_p_sums{$analysis}->{$family}->{$partition}
                  ->{S_diffs} += $line_arr[ $header_indices{'S_diffs'} ];
            }
        }
        close IN_BETWEEN_GROUP_CODON_FILE;

        my $between_product_dS;
        if ( $between_product_S_sites_sum > 0 ) {
            $between_product_dS =
              $between_product_S_diffs_sum / $between_product_S_sites_sum;
        }
        else {
            $between_product_dS = '*';
        }

        ##################################################################################
 # LOOP THROUGH EACH ANALYSIS, FAMILY, AND PARTITION FOR SLIDING WINDOW ANALYSIS
        foreach
          my $this_analysis ( sort keys %between_group_partition_codon_data_hh )
        {
            foreach my $this_family ( sort { $a <=> $b }
                keys %{ $between_group_partition_codon_data_hh{$this_analysis} }
              )
            {
                foreach my $this_partition (
                    sort { $a <=> $b } keys %{
                        $between_group_partition_codon_data_hh{$this_analysis}
                          ->{$this_family}
                    }
                  )
                {

                    # Get group names
                    my $group_name_i =
                      $between_group_a_f_p_groups{$this_analysis}
                      ->{$this_family}->{$this_partition}->{group_1};
                    my $group_name_j =
                      $between_group_a_f_p_groups{$this_analysis}
                      ->{$this_family}->{$this_partition}->{group_2};

#					print "\nTEST RETRIEVAL:\nanalysis: $this_analysis\nfamily: $this_family\npartition: $this_partition\n\n";

                    #############################################################
                    # SLIDING WINDOW WITH BOOTSTRAPS
                    # Generate and store sliding window and bootstrap analyses
                    my %sliding_window_results_hh
                      ; # Will store first_codon->N_diffs_window/S_diffs_window/...

#	my $sw_file_name = "partition$partition_id\_sw_" . $sliding_window_size . "codons\_results.txt";
#
#	if (-e "$sw_file_name") {
#		die "\n\n### WARNING:\n## The file $sw_file_name already exists.\n## Please ".
#			"rename or move this file so a new one ".
#			"can be created.\n\n";
#	}

                    my @codons = keys(
                        %{
                            $between_group_partition_codon_data_hh{
                                $this_analysis}->{$this_family}
                              ->{$this_partition}
                        }
                    );
                    my @sorted_codons     = sort { $a <=> $b } (@codons);
                    my $num_sorted_codons = scalar @sorted_codons;

                    #					print "\n\ncodons: @sorted_codons\n\n";

                    if ( $sliding_window_size > $num_sorted_codons ) {
                        die
"\n\n### WARNING:\n## For partition $this_partition\,\n## there are "
                          . "not enough codons to perform a sliding window of $sliding_window_size.\n\n";
                    }

                    # FOR EACH WINDOW
                    my $window_num  = 0;
                    my $first_print = 0;

                    #					print "\nProcessing sliding windows...\n";

# CODON_INDEX must begin one before the smallest actual codon number... don't assume it's 1
                    for (
                        my $codon_index = 0 ;
                        $codon_index <=
                        ( $num_sorted_codons - $sliding_window_size ) ;
                        $codon_index += $sliding_window_step
                      )
                    {
                        $window_num++;

                        my $first_codon = $codon_index + $sorted_codons[0];
                        my $last_codon =
                          $first_codon + $sliding_window_size - 1;

#						print "#### PROCESSING WINDOW $window_num\, codons $first_codon\-$last_codon\... \n";

                        if ( $last_codon <=
                            ( $num_sorted_codons + $sorted_codons[0] ) )
                        {    # as long as we're in range

                            my $window_N_diffs;
                            my $window_S_diffs;
                            my $window_N_sites;
                            my $window_S_sites;

                            for (
                                my $i = $first_codon ;
                                $i <= $last_codon ;
                                $i++
                              )
                            {
#								print "\n\nTo codon $i, we add: " . $between_group_partition_codon_data_hh{$this_analysis}->{$this_family}->{$this_partition}->{$i}->{N_diffs} . "\n\n";
#								print "\nTEST RETRIEVAL:\nanalysis: $this_analysis\nfamily: $this_family\npartition: $this_partition\n\n";
                                $window_N_diffs +=
                                  $between_group_partition_codon_data_hh{
                                    $this_analysis}->{$this_family}
                                  ->{$this_partition}->{$i}->{N_diffs};
                                $window_S_diffs +=
                                  $between_group_partition_codon_data_hh{
                                    $this_analysis}->{$this_family}
                                  ->{$this_partition}->{$i}->{S_diffs};
                                $window_N_sites +=
                                  $between_group_partition_codon_data_hh{
                                    $this_analysis}->{$this_family}
                                  ->{$this_partition}->{$i}->{N_sites};
                                $window_S_sites +=
                                  $between_group_partition_codon_data_hh{
                                    $this_analysis}->{$this_family}
                                  ->{$this_partition}->{$i}->{S_sites};
                            }

                            my $window_dN;
                            if ( $window_N_sites > 0 ) {
                                $window_dN =
                                  ( $window_N_diffs / $window_N_sites );
                            }
                            else {
                                $window_dN = '*';
                            }

                            my $window_dS;
                            if ( $window_S_sites > 0 ) {
                                $window_dS =
                                  ( $window_S_diffs / $window_S_sites );
                            }
                            else {
                                $window_dS = '*';
                            }

                            my $window_w;
                            if ( $window_dS > 0 && $window_dS ne '*' ) {
                                $window_w = ( $window_dN / $window_dS );
                            }
                            else {
                                $window_w = '*';
                            }

                            my $window_dN_minus_dS;
                            if ( $window_dN >= 0 && $window_dS >= 0 ) {
                                $window_dN_minus_dS = $window_dN - $window_dS;
                            }

                            #		print "\n\n";

                            #################################################################################
                            # BOOTSTRAP the WINDOW
                            my $SE_dN_minus_dS = 0;
                            my $window_boot_Z  = 'NA';

                            if ( $num_bootstraps > 1 ) {

                                # MAKE BOOTSTRAP FILE

                       #print "num_codons_in_product: $num_codons_in_product\n";

                                my @sim_N_sites_arr;
                                my @sim_S_sites_arr;
                                my @sim_N_diffs_arr;
                                my @sim_S_diffs_arr;

                                for (
                                    my $bootstrap_num = 1 ;
                                    $bootstrap_num <= $num_bootstraps ;
                                    $bootstrap_num++
                                  )
                                {
                                    srand();

  # New bootstrap run, of $num_bootstraps
  #if($this_family == 460) {
  #	print "\nPartition $this_partition \| BOOTSTRAP REPLICATE $bootstrap_num\n";
  #}

                                    my $sim_N_sites_sum = 0;
                                    my $sim_S_sites_sum = 0;
                                    my $sim_N_diffs_sum = 0;
                                    my $sim_S_diffs_sum = 0;

# METHOD
# (1) loop through $num_codons_in_product
# (2) for each, select a random codon position's results: Nd, Sd, Ns, Ss.
# where are the aforementioned stored? As before:
# $between_group_partition_codon_data_hh{$this_analysis}->{$this_family}->{$this_partition}->{$codon_num}->{N_diffs}, etc.

               # SAMPLE codon sites, up to the actual number of codons in window
                                    for (
                                        my $window_codon_index = 0 ;
                                        $window_codon_index <
                                        $sliding_window_size ;
                                        $window_codon_index++
                                      )
                                    {
        # Choose a random codon between this window $first_codon and $last_codon
                                        my $random_codon_num = int(
                                            rand(
                                                $last_codon - $first_codon + 1
                                            )
                                          ) +
                                          $first_codon;

                                        my $sampled_N_diffs =
                                          $between_group_partition_codon_data_hh{
                                            $this_analysis}->{$this_family}
                                          ->{$this_partition}
                                          ->{$random_codon_num}->{N_diffs};
                                        my $sampled_S_diffs =
                                          $between_group_partition_codon_data_hh{
                                            $this_analysis}->{$this_family}
                                          ->{$this_partition}
                                          ->{$random_codon_num}->{S_diffs};
                                        my $sampled_N_sites =
                                          $between_group_partition_codon_data_hh{
                                            $this_analysis}->{$this_family}
                                          ->{$this_partition}
                                          ->{$random_codon_num}->{N_sites};
                                        my $sampled_S_sites =
                                          $between_group_partition_codon_data_hh{
                                            $this_analysis}->{$this_family}
                                          ->{$this_partition}
                                          ->{$random_codon_num}->{S_sites};

                                        $sim_N_diffs_sum += $sampled_N_diffs;
                                        $sim_S_diffs_sum += $sampled_S_diffs;
                                        $sim_N_sites_sum += $sampled_N_sites;
                                        $sim_S_sites_sum += $sampled_S_sites;

                                    } # finished compiling all sampled codons (cols in alignment)

                                    push( @sim_N_diffs_arr, $sim_N_diffs_sum );
                                    push( @sim_S_diffs_arr, $sim_S_diffs_sum );
                                    push( @sim_N_sites_arr, $sim_N_sites_sum );
                                    push( @sim_S_sites_arr, $sim_S_sites_sum );

                                    #print $random_sum . "\n";

                                }    # end last bootstrap

                                # Recall we have $window_dN_minus_dS above

            # CALCULATE BOOTSTRAP STANDARD ERROR HERE; NEI & KUMAR (2000) METHOD
                                my @sim_dN_minus_dS;

                                for (
                                    my $sim_num = 0 ;
                                    $sim_num < scalar(@sim_N_sites_arr) ;
                                    $sim_num++
                                  )
                                {
                                    my $this_round_dN = '*';
                                    if ( $sim_N_sites_arr[$sim_num] > 0 ) {
                                        $this_round_dN =
                                          $sim_N_diffs_arr[$sim_num] /
                                          $sim_N_sites_arr[$sim_num];
                                    }

                                    my $this_round_dS = '*';
                                    if ( $sim_S_sites_arr[$sim_num] > 0 ) {
                                        $this_round_dS =
                                          $sim_S_diffs_arr[$sim_num] /
                                          $sim_S_sites_arr[$sim_num];
                                    }

                                    my $this_round_dN_minus_dS = '*';
                                    if (   $this_round_dN >= 0
                                        && $this_round_dS >= 0 )
                                    {
                                        $this_round_dN_minus_dS =
                                          $this_round_dN - $this_round_dS;
                                    }

                                    push( @sim_dN_minus_dS,
                                        $this_round_dN_minus_dS );
                                }

                                $SE_dN_minus_dS =
                                  &standard_deviation(@sim_dN_minus_dS);

                                if ( $SE_dN_minus_dS > 0 ) {
                                    $window_boot_Z =
                                      $window_dN_minus_dS / $SE_dN_minus_dS;
                                }

                            }    # END BOOTSTRAPS for WINDOW

                            my $partition_between_group_dS = '*';
                            if ( $between_group_a_f_p_sums{$this_analysis}
                                ->{$this_family}->{$this_partition}->{S_sites}
                                > 0 )
                            {
                                $partition_between_group_dS =
                                  $between_group_a_f_p_sums{$this_analysis}
                                  ->{$this_family}->{$this_partition}->{S_diffs}
                                  / $between_group_a_f_p_sums{$this_analysis}
                                  ->{$this_family}->{$this_partition}
                                  ->{S_sites};
                            }

                            my $window_dN_exceeds_product_dS = '';
                            if ( $window_dN > $between_product_dS ) {
                                $window_dN_exceeds_product_dS = '*';
                            }

                            my $min_defined_codons_in_window =
                              1000;    # or something big

                            my $window_num_defined_codons_g1 = 0;
                            for (
                                my $j = $first_codon ;
                                $j <= $last_codon ;
                                $j++
                              )
                            {
                                $window_num_defined_codons_g1 +=
                                  $between_group_partition_codon_data_hh{
                                    $this_analysis}->{$this_family}
                                  ->{$this_partition}->{$j}
                                  ->{'num_defined_codons_g1'};

                                if (
                                    $between_group_partition_codon_data_hh{
                                        $this_analysis}->{$this_family}
                                    ->{$this_partition}->{$j}
                                    ->{'num_defined_codons_g1'} <
                                    $min_defined_codons_in_window
                                  )
                                {
                                    $min_defined_codons_in_window =
                                      $between_group_partition_codon_data_hh{
                                        $this_analysis}->{$this_family}
                                      ->{$this_partition}->{$j}
                                      ->{'num_defined_codons_g1'};
                                }
                            }

                            my $window_num_defined_codons_g2 = 0;
                            for (
                                my $j = $first_codon ;
                                $j <= $last_codon ;
                                $j++
                              )
                            {
                                $window_num_defined_codons_g2 +=
                                  $between_group_partition_codon_data_hh{
                                    $this_analysis}->{$this_family}
                                  ->{$this_partition}->{$j}
                                  ->{'num_defined_codons_g2'};

                                if (
                                    $between_group_partition_codon_data_hh{
                                        $this_analysis}->{$this_family}
                                    ->{$this_partition}->{$j}
                                    ->{'num_defined_codons_g2'} <
                                    $min_defined_codons_in_window
                                  )
                                {
                                    $min_defined_codons_in_window =
                                      $between_group_partition_codon_data_hh{
                                        $this_analysis}->{$this_family}
                                      ->{$this_partition}->{$j}
                                      ->{'num_defined_codons_g2'};
                                }
                            }

                            if (   $min_defined_codons_in_window == 0
                                || $min_defined_codons_in_window eq ''
                                || !defined($min_defined_codons_in_window) )
                            {
                                $min_defined_codons_in_window = 0;
                            }

                            my $window_mean_num_defined_codons =
                              ( $window_num_defined_codons_g1 +
                                  $window_num_defined_codons_g2 ) / 2;

                            my $out_line_between_sw =
"$this_analysis\t$this_family\t$this_partition\t$group_name_i\t$group_name_j\t"
                              . "$window_num\t$first_codon\t$last_codon\t"
                              . "$window_num_defined_codons_g1\t$window_num_defined_codons_g2\t$window_mean_num_defined_codons\t$min_defined_codons_in_window\t"
                              . "$window_N_sites\t$window_S_sites\t$window_N_diffs\t$window_S_diffs\t"
                              . "$window_dN\t$window_dS\t$window_dN_minus_dS\t$window_w\t$window_dN_exceeds_product_dS";

                            if ( $num_bootstraps > 1 ) {

                                my $significance = '';

                                if ( $SE_dN_minus_dS > 0 ) {

                                    if ( abs($window_boot_Z) > 2.81 ) {
                                        $significance = '***';
                                    }
                                    elsif ( abs($window_boot_Z) > 1.96 ) {
                                        $significance = '**';
                                    }
                                    elsif ( abs($window_boot_Z) > 1.64 ) {
                                        $significance = '*';
                                    }

                                }

                                $out_line_between_sw .=
"\t$SE_dN_minus_dS\t$window_boot_Z\t$significance\n"
                                  ;    #####

                            }
                            else {
                                $out_line_between_sw .= "\n";
                            }

                            open( OUTFILE_BETWEEN_SW,
                                    ">>between\_group\_sw_"
                                  . $sliding_window_size
                                  . "codons\_results.txt" );
                            print OUTFILE_BETWEEN_SW "$out_line_between_sw";
                            close OUTFILE_BETWEEN_SW;

                        }
                    }
                }    # done looping partitions in the family
            }    # done looping families in the analysis
        }    # done looping analyses in the file
    }    # finish case in which between_group_codon_results.txt exists

    $pm_family->finish;    # PARALLELIZATION

}    # done looping FAMILIES

$pm_family->wait_all_children;    # PARALLELIZATION

print "\n\n";

#########################################################################################
# END PROGRAM and print a completion message to screen
&end_the_program;

#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################

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
sub get_header_names {

# Originally assumed that we've received a tempfile ending in "_snpg9temp.txt"
# However, we're now calling it at least once before creating the tempfile to
# see what kind of processing (e.g., Geneious to CLC) is needed prior to tempfile
# creation. Must include capability to get headers for .CSV file
    my ( $curr_snp_report_filename, $filename ) = @_;

    #print "\n$curr_snp_report_filename\n";

    #my $newline_char = &detect_newline_char($curr_snp_report_filename);
    #my $old_newline = $/;
    #$/ = $newline_char;

    my $seen_tab_delimited     = 0;
    my $seen_comma_delimited   = 0;
    my $seen_vcf_tab_delimited = 0;
    my @line_arr;

    my $line = 0;
    open( CURRINFILE, $curr_snp_report_filename );

    #seek(CURRINFILE,0,0);
    while (<CURRINFILE>) {

        #print "$_";
        if ( $line == 0 ) {
            chomp;

            # CHOMP for 3 operating systems
            if ( $_ =~ /\r\n$/ ) {
                $_ =~ s/\r\n//;
            }
            elsif ( $_ =~ /\r$/ ) {
                $_ =~ s/\r//;
            }
            elsif ( $_ =~ /\n$/ ) {
                $_ =~ s/\n//;
            }

            if ( $_ =~ /\t\w+\t/ ) {    # it's TAB-delimited
                @line_arr = split( /\t/, $_ );

                #print "TAB!!!!!";
                last;
            }
            elsif ( $_ =~ /,\w+,/ ) {    # it's COMMA-delimited
                @line_arr = split( /,/, $_ );

                #print "COMMA!!!!!";
                last;
            }

            $line++;
        }
        elsif ( $line > 0 && $_ =~ /^##/ ) {
            $line++;
        }
        elsif ( $line > 0 && ( $_ =~ /^#CHROM/ ) ) {
            chomp;

            # CHOMP for 3 operating systems
            if ( $_ =~ /\r\n$/ ) {
                $_ =~ s/\r\n//;
            }
            elsif ( $_ =~ /\r$/ ) {
                $_ =~ s/\r//;
            }
            elsif ( $_ =~ /\n$/ ) {
                $_ =~ s/\n//;
            }

            if ( $_ =~ /\t/ ) {    # it's TAB-delimited
                @line_arr = split( /\t/, $_ );

                #print "TAB!!!!!";
                last;
            }
            elsif ( $_ =~ /,/ ) {    # it's COMMA-delimited
                @line_arr = split( /,/, $_ );

                #print "COMMA!!!!!";
                last;
            }
            else {
                chdir('SNPGenie_Results');
                open( ERROR_FILE, ">>SNPGenie\_LOG\.txt" );

                # FILE | PRODUCT | SITE | CODON | WARNING

            # No change OR error should occur if the file does not, in fact, end
            # with this SUFFIX
                my $file_nm = $curr_snp_report_filename;

                #$file_nm =~ s/_snpg9temp.txt/.txt/;
                $file_nm =~ s/_\w\w\w\w.txt/.txt/;

                print ERROR_FILE "$filename\tN/A\tN/A\t"
                  . "File not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPGenie terminated.\n";
                close ERROR_FILE;
                chdir('..');

                #unlink $curr_snp_report_filename;

                die "\n\n## WARNING: The SNP Report $filename is "
                  . "not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPgenie "
                  . "terminated\n\n";
            }
        }
        else {
            chdir('SNPGenie_Results');
            open( ERROR_FILE, ">>SNPGenie\_LOG\.txt" );

            # FILE | PRODUCT | SITE | CODON | WARNING

            # No change OR error should occur if the file does not, in fact, end
            # with this SUFFIX
            my $file_nm = $curr_snp_report_filename;

            #$file_nm =~ s/_snpg9temp.txt/.txt/;
            $file_nm =~ s/_\w\w\w\w.txt/.txt/;

            print ERROR_FILE "$filename\tN/A\tN/A\t"
              . "File not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPGenie terminated.\n";
            close ERROR_FILE;
            chdir('..');

            #unlink $curr_snp_report_filename;

            die "\n\n## WARNING: The SNP Report $filename is "
              . "not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPgenie "
              . "terminated\n\n";
        }
    }
    seek( CURRINFILE, 0, 0 );
    close CURRINFILE;

    #$/ = $old_newline;
    return @line_arr;
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
      . "\n##               SNPGenie for OrthologID completed successfully.              ##"
      . "\n##        Please find results in the \/data\/ dir and subdirs (families).       ##\n"
      . "################################################################################"
      . "\n\n\n";
}
