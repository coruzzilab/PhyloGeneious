#!/usr/bin/env perl

use strict;
use Cwd;
use Time::Piece;
my @pidtbl  = ();
my @opentbl = ();

#my @tntcmd = qw(oid2.proc oid3.proc);
my @tntcmd = ();
my $nproc  = 2;
my $tntstr = 'pf';
my $nthr   = 2;
my $curthr = 0;
my $tm;
my $savDir = getcwd;
$tm = localtime;
print "runtnt in $savDir at $tm\n";
print "runtntmx in $savDir\n";

#    chdir "1";
die "runtntmx.pl nproc nthreads prefix\n" if ( @ARGV < 3 );
$nproc  = $ARGV[0];
$nproc  = 2 if $nproc < 2;
$nthr   = $ARGV[1];
$nthr   = $nproc if ( $nthr > $nproc );
$tntstr = $ARGV[2];
print "run $nproc procs max threads $nthr, prefix= $tntstr\n";

for ( my $i = 1 ; $i <= $nproc ; $i++ ) {
    my $procn = $i + 1;
    push @tntcmd, "$tntstr$procn.proc";
}

#    print join ",",@tntcmd,"\n";

foreach my $tnt (@tntcmd) {
    my $pid;
    sleep 2;

    #        next if $pid = fork;  #parent
    if ( $pid = fork() ) {    # parent
        $curthr++;
        next if ( $curthr < $nthr );
        my $child = wait();
        $tm = localtime;
        print "$tm $child completed\n";
        $curthr--;
        next;    #so we can fork
    }

    $tm = localtime;
    print "$tm starting $tnt\n";
    my $status = system("tnt p $tnt 0</dev/null 1>&0 2>&0");
    print "tnt with $tnt exit $status\n";
    exit;
}
1 while ( wait() != -1 );

#    my $ov;
#    $pid = open($ov,"tnt.command p tr1.proc 0</dev/null |") or print "first fork fails $1\n";
#    push @pidtbl,$pid;
#    push @opentbl,$ov;
#    $pid = open($ov,"tnt.command p tr2.proc  |") or print "second fork fails $1\n";
#    push @pidtbl,$pid;
#    push @opentbl,$ov;
#    foreach $pid(@pidtbl){
#        my $vol = waitpid($pid,0);
#        print "$pid return $vol\n";
#    }
#    foreach $ov(@opentbl) {
#        close $ov;
#    }
my $one = '1';
my $tnt = "$tntstr$one.proc";
$tm = localtime;
print "$tm starting $tnt\n";
my $status = system("tnt p $tnt 0</dev/null 1>&0 2>&0");
$tm = localtime;
print "$tm tnt with $tnt exit $status\n";

my $lsv = `ls -l $tntstr*.tre`;
print "ls value $lsv\n";

