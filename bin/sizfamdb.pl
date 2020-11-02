#!/usr/bin/env perl
#
# OrthologID main script
#
#

use warnings;

my $OID_HOME;
my $OID_USER_DIR;
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
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
use DB_File;
use Time::Piece;
use strict;
my $OID_DATADIR    = "$OID_USER_DIR/data";
my $familyDB = "$OID_USER_DIR/familyDB.dbm";
my %sztaxa = ();  #this is a dictionary of the size of all families in current data directory
my $cntren = 0;
tie (%sztaxa,"DB_File",$familyDB,O_WRONLY|O_CREAT,0666) or
   die "could not create tie $familyDB\n";

	    chdir $OID_DATADIR;
	foreach my $dir (<[1-9]*>) {
        chdir $dir;
        my $pid = open(GR, "grep -c '>' FAMILY |") ;
        $_ = <GR>;
        chomp;
        my $dirsize = $_;
        $sztaxa{$dir} = $_ if $dirsize>3;  #only use families size>3
	    chdir $OID_DATADIR;
        if ($dirsize < 4) {
            my $newdir = "S"."$dir";
            rename ($dir,$newdir);
            $cntren++;
        }
    }
    my @kefnd = keys %sztaxa;
    my $fam = scalar @kefnd;
    print "there are $fam families renamed $cntren small families\n";
