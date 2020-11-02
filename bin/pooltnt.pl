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
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
    $HPC=$ENV{'HPC'};
    $HPC='P' if ! defined($HPC);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);
    $OID_BLAST = "$OID_USER_DIR/blast";
}

use lib "$OID_HOME/lib";
use OrthologID; 
my @pidtbl = ();
my @opentbl = ();
#my @tntcmd = qw(oid2.proc oid3.proc);
my @tntcmd = ();
my $nproc = 2;
my $tntstr = 'pf';
my $nthr = 2;
my $curthr = 0;
my $tm;
my $pooldone = 'log/job/pooldone';
	my $savDir = getcwd;
    $tm = localtime;
    print "pooltnt.pl in $savDir at $tm\n";
#    chdir "1";
#    die "runtntmx.pl nproc nthreads prefix\n" if (@ARGV<3);
    die "pooltnt.pl startfam,endfam,ncpu\n" if (@ARGV<3);
    my $startfam = $ARGV[0];
    my $endfam = $ARGV[1];
        $nthr = $ARGV[2];
        print "we run from fam $startfam to $endfam on $nthr cpus\n";


#    foreach my $tnt (@tntcmd){
    for (my $fam=$startfam; $fam <=$endfam;$fam++){ 
        my $pid;
#        sleep 2;
#        next if $pid = fork;  #parent
        if ($pid = fork()) { # parent
            $curthr ++;
            next if ($curthr < $nthr);
            my $child = wait();
            $tm = localtime;
            print "$tm $child completed\n";
            $curthr--;
            next; #so we can fork
        }


            $tm = localtime;
        print "$tm starting $fam\n";
        if (-s "data/$fam/oid.tre"){
             print "$fam already done\n";
             exit;
         }
         my $famdol = '^'."$fam".'$';
         alignFamily("$famdol");
         makeTree("$famdol");
#        my $status = system("tnt p $tnt 0</dev/null 1>&0 2>&0");
#
        my $tm = localtime;
        print "$tm tnt with $fam done \n";
        exit;
    }
    1 while (wait() != -1);

    open (FI,">$pooldone");
    print FI "$startfam,$endfam\n";
    close FI;
    $tm = localtime;
    print "$tm pooltnt done\n";
    
