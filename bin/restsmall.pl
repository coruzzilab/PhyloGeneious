#!/usr/bin/env perl
#
# OrthologID main script
#
#

use warnings;

my $OID_HOME;
my $OID_USER_DIR;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
}

use lib "$OID_HOME/lib";
use OrthologID;
use Getopt::Std;
use Cwd;
use Fcntl;
use AnyDBM_File;
use Time::Piece;
use strict;
my $OID_DATADIR = "$OID_USER_DIR/data";
my $cntren      = 0;
my $cntfnd      = 0;

chdir $OID_DATADIR;
foreach my $dir (<S[1-9]*>) {
    chdir $dir;
    my $pid = open( GR, "grep -c '>' FAMILY |" );
    $_ = <GR>;
    chomp;
    $cntfnd++;
    my $dirsize = $_;
    chdir $OID_DATADIR;
    if ( $dirsize < 4 ) {
        my $newdir = substr( $dir, 1 );    #skip "S"
        rename( $dir, $newdir );
        $cntren++;
    }
}
print "renamed $cntren small families ($cntfnd found)\n";
