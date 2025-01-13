#!/usr/bin/env perl

use strict;
use Cwd;
use Time::Piece;

use Fcntl;
my $OID_HOME;
my $OID_USER_DIR;
my $OID_BLAST;
my $OID_BLASTDIR;
my $HPC;
my $fam;
my $tm;
my $procstat;
#my $alignedFamily;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $HPC          = $ENV{'HPC'};
    $HPC          = 'P' if !defined($HPC);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
    $OID_BLAST = "$OID_USER_DIR/blast";
}

use lib "$OID_HOME/lib";
use OrthologID;

$fam  = $ARGV[0];

$tm = localtime;
print "$tm starting $fam\n";
if ( -s "data/$fam/oid.tre" ) {
    print "$fam already done\n";
    exit;
}

$procstat = system("grep -i -q 'DNA' $OID_USER_DIR/procfiles.txt");
if ( $procstat > 0 ) { #re-set procfiles
    system("cp $OID_HOME/testdata/procfiles_dna.txt $OID_USER_DIR/procfiles.txt");
    OrthologID::init_proc();
}

$OrthologID::alignedFamily = "FAMILY.aligned.revfasta";
local $SIG{ALRM} = sub { print "timed out $fam\n"; exit };
alarm 600;    #these are fast rtns but may hang
makeTree("$fam");
alarm 0;      #cancel alarm it works

my $tm = localtime;
print "$tm tnt with $fam done \n";
exit;
