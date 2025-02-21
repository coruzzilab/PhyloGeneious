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
use DB_File;
use Time::Piece;
use strict;
my $OID_DATADIR = "$OID_USER_DIR/data";
my $familyDB    = "$OID_USER_DIR/familyDB.dbm";
my %sztaxa      = ()
  ;  #this is a dictionary of the size of all families in current data directory

if ( !-s "$familyDB" ) {
    if (! $ENV_WRAPPER eq "") {
        my $rc = system("$OID_HOME/bin/sizfamdb.pl"); #$ENV_WRAPPER 
    }
    else {
        my $rc = system("$OID_HOME/bin/sizfamdb.pl");
    }
}
my $badtie  = 0;
my $cntmiss = 0;
tie( %sztaxa, "DB_File", $familyDB, O_RDONLY, 0666 )
  or do { print "could not create tie $familyDB\n"; $badtie = 1; };
my @taxsz = keys %sztaxa;
my $ntaxa = @taxsz;
print "size of taxa table  $ntaxa\n";
chdir $OID_DATADIR;
my $fndtx = 0;

if ( !$badtie ) {
    foreach my $dir (<[1-9]*>) {
        if ( $sztaxa{$dir} ) {
            $fndtx++;
        }
        else {
            print "data/$dir not in sztaxa\n" if $cntmiss < 20;
            $cntmiss++;
        }
    }

    untie %sztaxa;
}    #if any taxa
BADTIE:
if ( $cntmiss > 0 or $badtie ) {
    print "somthing wrong with $familyDB\n";
    if ( exists $ARGV[0] ) {
        print "already retried - exit \n";
        exit(1);
    }
    unlink "$familyDB";
    unlink "$familyDB.db";    #on mac
    if (! $ENV_WRAPPER eq "") {
        my $rc = system("$OID_HOME/bin/rdfamdb.pl retry"); #$ENV_WRAPPER 
        exit($rc);
    }
    else {
        my $rc = system("$OID_HOME/bin/rdfamdb.pl retry");
        exit($rc);
    }
}
exit(0);

