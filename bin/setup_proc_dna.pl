#!/usr/bin/env perl
## this is a paraphrase ofdameon's script to do a lot of rathchets
use Time::Piece;
use Time::Seconds;
use Cwd;
use Getopt::Long;
use POSIX qw/ceil/;

my $OID_HOME;
my $OID_USER_DIR;
my $OID_BIN;
my $ENV_WRAPPER;
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);
    $OID_BIN="$OID_HOME/bin";
    $ENV_WRAPPER  = $ENV{'ENV_WRAPPER'};
    if (!defined($ENV_WRAPPER)){
        $ENV_WRAPPER = "";
    }
}

use feature "switch";
use strict;
my $t = localtime;
my $fam;
#my $nchar;
#my $ntax;

    die "setup_proc_dna.pl family\n" if (@ARGV<1);
    $fam = $ARGV[0];
#    $ntax = $ARGV[1];
#    $nchar = $ARGV[2];

my $rtnum = 2;
my $time= "00:30:00";
my $prefix = "oi";
my $ntsk = 50;
my $nthr = 12;
my $iter;
#high rat only
    $rtnum = 6;
    $iter = 1;
    $prefix = "pi";
    $ntsk = 8;
    $nthr = 4;
    my $rc;
my $rtn = "$OID_BIN/ratbld"."$rtnum"."_dna.pl";
my $timeout;
    if (! $ENV_WRAPPER eq "") {
        my $rc = system("perl $rtn -n $ntsk  -m FAMILY.aligned.revfasta -i $iter -h 10 -p $prefix -t 4"); #$ENV_WRAPPER 
    }
    else {
        $rc = system("perl $rtn -n $ntsk  -m FAMILY.aligned.revfasta -i $iter -h 10 -p $prefix -t 4");
    }
    $timeout = "96:00:00";
open(FB,">setup.txt");
print FB  "$ntsk\t$nthr\t$prefix\t$timeout\n";
close FB
