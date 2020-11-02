#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;
use Fcntl;
my $OID_HOME;
my $OID_USER_DIR;
my $OID_DATADIR ;
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);
    $OID_DATADIR = "$OID_USER_DIR/data";
}

use lib "$OID_HOME/lib";
use OrthologID; 

sub makenex {

	my $dirRE = shift;
    if (! defined ($dirRE)){
        print "aligned family called withot argument\n";
        return;
    }
#	$dirRE = '.*' if !defined($dirRE);
    my $dir = $dirRE;
    if ($dirRE =~ /\^(\d+)\$/){
        $dir = $1;
    }
	my $verbose = $OrthologID::VERBOSE;
    my $alignedFamily = $OrthologID::alignedFamily;
    my $verbose = $OrthologID::verbose;
    print "called makenex.pl $dir $alignedFamily\n";

	# Defaults
	my $prefix    = "oid";
	#my $paup_exe  = "paup";
	my $tnt_exe = "tnt";
	my $nrat      = 10;       # number of ratchets
	my $niter     = 100;      # number of iterations per ratchet
	my $timeLimit = 10800;    # Time limit for a heuristic search
	my $ratTimeLimit = 60;    # Time limit for a search within a ratchet iteration

	my $thisFunc = ( caller(0) )[3];
#	my $thisFunc = 'makenex';
	my $oldFH;
	my %seq;
	my @outgroup;           # Outgroup taxa for current family
    my $curcut = 0;
    my $maxcut = 0;
    my $curproc = 0;
    my $ntaxa = 0;


    $verbose = 2;
	# Change to data dir
	my $savDir = getcwd;
	chdir $OID_DATADIR;

	# Go over each family
	my $count = 0;
	print "Making trees ...\n" if $verbose;
#	foreach my $dir (<[1-9]*>) 
#		next if $dir !~ /$dirRE/;
    for (my $x=0;$x<1;$x++){
		$count++;
		chdir $dir
		  or die "$thisFunc: failed to change to family directory $dir: $!\n";
		if ( $verbose == 1 ) {
			$| = 1;    # flush STDOUT
			print "\b" x 7;
			printf "%6d", $count;
			$| = 0;
		}
		elsif ( $verbose == 2 ) {
			print "Tree search for family $dir\n";
		}

		%seq      = ();
		@outgroup = ();

		# Check nexus file or tree file does not already exist
		if ( -s $prefix . ".tre" ) {
			print "$dir Tree file already exists ... skipping\n"
			  if $verbose;
#			chdir $OID_DATADIR;
#			next;
	        chdir $savDir;
            return;
		}

		open( ALIGN, "<", $alignedFamily )
		  or die "$thisFunc: failed to open alignment: $!\n";
		my $nam;
		while (<ALIGN>) {
			chomp;
			if (/^>/) {
				$nam = substr( $_, 1 );
				$seq{$nam} = "";
			}
			elsif ( defined $seq{$nam} ) {
				$seq{$nam} .= $_;
			}
		}
		close ALIGN;
		my $length = length( $seq{$nam} ); # assume they are all the same length
		my $ntaxa  = keys %seq;
		if ($ntaxa <4){
			print "$dir Too few taxa. TNT cannot create tree with less than 4 taxa\n"
			 if $verbose;
#			chdir $OID_DATADIR;
#			next;
	        chdir $savDir;
            return;
		}
#		if (-s $prefix.".nex") goto RUNTNT 
		open( NEX, ">", $prefix . ".nex" )
		  or die "$thisFunc: Failed to open file for writing: $!\n";
		$oldFH = select(NEX);
		print "xread\n";
		print "$length $ntaxa\n";
		print "&[proteins nogaps]\n";
		foreach ( sort keys %seq ) {
			if (length($_) < 32) {
				# Try to line up the characters
				print $_ . " " x ( 32 - length($_) );
			}
			else {
				print "$_ ";
			}
			print $seq{$_} . "\n";
			# Check if taxon is outgroup
			foreach my $og (@OrthologID::OUTGROUP) {
				if (/^$og#/) {
					push @outgroup, $_;
					last;
				}
			}
		}
		print ";\n";
		close(NEX);
		select($oldFH);
    }
}
#### dummy main
    die "makenex family \n" if (@ARGV<1);
    my $fam = $ARGV[0];
    makenex($fam);
