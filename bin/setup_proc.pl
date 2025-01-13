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
my $OID_WRAPPER;
my $TREEPROGRAM;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
    $OID_WRAPPER = $ENV{'OID_WRAPPER'};
    $OID_BIN = "$OID_HOME/bin";
}

use feature "switch";
use strict;
my $t = localtime;
my $fam;
my $nchar;
my $ntax;
my $use_oblong=0;


die "setup_proc.pl family ntaxa nchar\n" if ( @ARGV < 3 );
$fam   = $ARGV[0];
$ntax  = $ARGV[1];
$nchar = $ARGV[2];
$TREEPROGRAM = $OrthologID::TREEPROGRAM;
if ($TREEPROGRAM =~ /OBLONG/) {
    $use_oblong = 1;
}

# for now very simplistic ntax < 1500 just run ratbld2.pl with 2 iter, < 2000 ratbld3.pl with 2 iter
# < 2500 ratbld3.pl with 1 iter and rest ratbld4.pl with only 1 iter
my $rtnum  = 2;
my $time   = "00:30:00";
my $prefix = "oi";
my $ntsk   = 50;
my $nthr   = 12;
my $iter;
if ( $ntax < 1500 ) {
    $rtnum = 2;
    if ( $nchar < 3000 ) {
        $iter = 2;
    }
    else {
        $iter = 1;
    }
    $prefix = "oi";
}
elsif ( $ntax < 2000 ) {
    $rtnum  = 3;
    $iter   = 3;
    $prefix = "op";
}
elsif ( $ntax < 2500 ) {
    $rtnum  = 3;
    $iter   = 2;
    $prefix = "op";
}
elsif ( $ntax < 5000 ) {
    $rtnum  = 3;
    $iter   = 1;
    $prefix = "op";
}
else {
    $rtnum  = 6;
    $iter   = 1;
    $prefix = "pi";
    $ntsk   = 8;
    $nthr   = 4;
}
my $rc;
my $rtn = "$OID_BIN/ratbld" . "$rtnum" . ".pl";
my $timeout;
if ($use_oblong == 1) {
    my $rtn = "$OID_BIN/oblong_ratbld2_subprocs.sbatch";
    $rc = system(
        "sbatch --time=$time -o logs/${prefix}.%J.log ${OID_WRAPPER} $rtn"
        );
    $timeout = "00:30:00";
}
else {
    if ( $rtnum < 4 ) {
    $rc = system(
        "$rtn -c $nchar -n $ntsk  -m oid.nex -i $iter -p $prefix -t 00:30:00");
    $timeout = "24:00:00";
    }
    elsif ( $rtnum < 6 ) {
        $rc = system("$rtn -n $ntsk  -m oid.nex -i $iter -h 10 -p $prefix -t 2");
        $timeout = "48:00:00";
    }
    else {
        $rc = system("$rtn -n $ntsk  -m oid.nex -i $iter -h 10 -p $prefix -t 4");
        $timeout = "96:00:00";
    }
}

open( FB, ">setup.txt" );
print FB "$ntsk\t$nthr\t$prefix\t$timeout\n";
close FB
