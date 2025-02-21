#!/usr/bin/env perl
#
# OrthologID main script
#
#

use warnings;

my $OID_HOME;
my $OID_USER_DIR;
my $ENV_WRAPPER;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
    $ENV_WRAPPER  = $ENV{'ENV_WRAPPER'};
    if (!defined($ENV_WRAPPER)){
        $ENV_WRAPPER = "";
    }
}

use lib "$OID_HOME/lib";
use OrthologID;
use Getopt::Std;
use Cwd;
use Fcntl;
use AnyDBM_File;
use Time::Piece;
use File::Copy;
use strict;
my $OID_DATADIR = "$OID_USER_DIR/data";
my $smalldir    = "$OID_DATADIR/smalldir";
my $cntren      = 0;
my $cntfnd      = 0;

if ( !-d $smalldir ) {
    mkdir $smalldir, 0755 or die "cannot make $smalldir";
}
chdir $OID_DATADIR;
foreach my $dir (<S[1-9]*>) {
    $cntfnd++;
    next if ( !-f "$dir/FAMILY" );
    my $sfilenm = 'F' . substr( $dir, 1 ) . '.faa';
    copy( "$dir/FAMILY", "$smalldir/$sfilenm" );
    my $alignam = substr( $sfilenm, 0, -3 ) . 'aligned';
    chdir $smalldir;
    if (! $ENV_WRAPPER eq "") {
        my $rc = system("$ENV_WRAPPER mafft --auto --quiet --anysymbol $sfilenm>$alignam");
    else {
        my $rc = system("mafft --auto --quiet --anysymbol $sfilenm>$alignam");
}

    chdir $OID_DATADIR;

    $cntren++;
}
print "created mafft align for $cntren small families ($cntfnd found)\n";
