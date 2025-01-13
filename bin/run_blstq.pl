#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;
use Fcntl;
my $OID_HOME;
my $OID_USER_DIR;
my $Blastdir;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
}

use lib "$OID_HOME/lib";
use Blstq;
use Getopt::Long;
use POSIX qw/ceil/;
use strict;
my $t     = localtime;
my @stuff = ();
my $strt;
my %qsid      = ();
my $blastpart = 0;
my $blstthr   = 0;
print "local time $t\n";

if ( @ARGV > 0 ) {
    $blastpart = "$ARGV[0]";
}
if ( @ARGV > 1 ) {
    $blstthr = "$ARGV[1]";
}
if ( $blastpart <= 0 ) {
    print "run_blastq.pl part [num_threads] \n";
    exit(1);
}
our $num_thr = $Blstq::num_thr;
$num_thr = $blstthr if ( $blstthr > 0 );
print "running blast part $blastpart with $num_thr threads\n";
runBlast($blastpart);
exit(0);

