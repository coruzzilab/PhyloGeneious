#!/usr/bin/env perl

use strict;
use Cwd;
my @pidtbl = ();
my @opentbl = ();
my @tntcmd = qw(tr1.proc tr2.proc);

	my $savDir = getcwd;
    chdir "1";
    foreach my $tnt (@tntcmd){
        my $pid;
        next if $pid = fork;  #parent
        my $status = system("tnt.command p $tnt 0</dev/null 1>&0 2>&0");
        print "tnt with $tnt exit $status\n";
        exit;
    }
    1 while (wait() != -1);

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
        $tnt = tr.proc;
        my $status = system("tnt.command p $tnt 0</dev/null 1>&0 2>&0");
        print "tnt with $tnt exit $status\n";
    
    my $lsv = `ls -l tr*`;
    print "ls value $lsv\n";



