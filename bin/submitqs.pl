#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;
use Fcntl;
my $OID_HOME;
my $OID_USER_DIR;
my $OID_BLAST;
my $NCPU;   # cpu's to give blast
my $gb;   #gb to run blast
my $MYUSER;
my $HPC;
my $HPCTY;
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);
    $HPC = $ENV{OID_HPC};
    $HPC = "PBS" if !defined($OID_HPC);
    if ($HPC =~ /PBS/) {
        $HPCTY=1;
    } elsif ($HPC =~ /SLU/) {
        $HPCTY=2
    } elsif ($HPC =~ /SGE/) {
        $HPCTY = 3
    }else {
        $HPCTY = 1;
        print "$HPC unknown using PBS\n";
        $HPC  = 'PBS';
    }
 
    $OID_BLAST="$OID_USER_DIR/blast";
    $MYUSER = `whoami`;
}

use lib "$OID_HOME/lib";
use OrthologID; 
use Getopt::Long;
use strict;

sub startqs { #starts a qsub
#                startqs(\@qsid,$blid,$wall);
#                startqs(\@qsid,\@req)
#                @req array has: mem,wall,nthr,JOB_SCRIPT,args
    (my $rqs, ,my $blid ) = @_;
    (my $mem, my $wall,my $nthr, my $JOB_SCRIPT) = $@blid[0..3];
    my $nargs = len($@blid)  - 4;
## now we format submit string based of HPC type
#
#    my $JOB_SCRIPT="$OID_USER_DIR/run_oid_job.sh";

    if (substr($JOB_SCRIPT,-2) == 'sh'){
        $JOB_SCRIPT = substr($JOB_SCRIPT,0,-2).'$HPC'.'.sh';
    }

    my $walmem = "mem=$gb".",walltime="."$wall";
    my $PARAMS = "-l nodes=1:ppn=$NCPU -l $walmem";
#    my $PARAMS = "-l nodes=1:ppn=$NCPU -l mem=8GB,walltime=12:00:00";
    my $submit = q/qsub  /."$PARAMS $JOB_SCRIPT".q/ -v arg1="-b"/.",arg2=$blid  2>/dev/null";
    print "$submit\n";
#    my $pid = open(QS,"qsub -l $PARAMS $JOB_SCRIPT -v arg1=/"$at/",arg2=/"^$fam.$dol/" 2>/dev/null |");
    my $pid = open(QS,"$submit |");
    my $t = localtime;
    my $qs = <QS>;
    if ($qs){
        chomp $qs;
        print "start qsid $qs for blastid $blid\n";
        $$rqs{$qs} = [$blid,$t,$wall]; #put $qs in its place
        return 1;
    } else {
        print "qsub somehow fails\n";
        ## bug in code - will never end
        return 0
    }
}  ## startqs

sub activqs { #test if qsid still running 
    my $rqs = shift;
    my %active = ();
    my @stuff = ();
    my $qs;
    my $fam;
    my %actfam = ();
    my $strt;
    my $wall;
    my $myjob = $ENV{"PBS_JOBID"};
    my $fndme = 0;
    my $flag = 1;
    $fndme = 1 if (! $myjob);   #non pbs system
    my $t = localtime;
#    print "testqs called @$rqs\n";
	$| = 1;    # flush STDOUT
#    print "enter testqs as $MYUSER\n";
    while ($flag) {
        $flag = 0;
    my $pid = open(QS,"qstat -e -u $MYUSER |");
    my $jnk = <QS>;
#    $jnk = <QS>;
#    print "testqs $jnk\n";
    while (<QS>){
        chomp;
        @stuff = split;
        next if ($stuff[0] =~ /\D/);
        print "qs $stuff[0], $stuff[-2]\n";
        if ($stuff[-2] =~ /R|Q/){
            if ($stuff[-2] =~ /Q/){
                print "$stuff[0] still on Q wait for start\n";
                $flag = 1;
                sleep 30;
                last;
            }
            $qs = "$stuff[0]";
            $fndme = 1 if ($qs == $myjob);
            next if ($qs == $myjob);
            my $hrs = substr ($stuff[-1],0,2);
            my $min = substr ($stuff[-1],3.2);
            my $sec = substr ($stuff[-1],5,2);
            $strt = $t - $hrs*ONE_HOUR - $min*ONE_MINUTE ;
            print "$stuff[-1] was ",$strt,"\n";
            
            $active{"$qs"} = [$strt, $stuff[-1],$stuff[-3]];  #put found pids in hash
#            print "$stuff[0] is active\n";
        }
    }
  }  ## while flag
    if (! $fndme) {
        print "did not find me job $myjob - something wrong - try later\n";
        return ($fndme,%actfam);
    }
    foreach $qs (keys %active){
        
#        my $pid = open(GR, "grep -c '>' FAMILY |") ;
#        $sztaxa{$dir} = <GR>;
        my $log = "log/job/*".$qs;
         print "$log\n";
         my $pid = open(GR, "grep 'BLASTing .. ' $log 2>/dev/null| ");
         $fam = 0;
         $_ = <GR>;
         next if (! $_);
         chomp;
         print "found $_\n";
         if (/\D+(\d+)/){
             $fam = $1;
         }
            print "active $qs found fam $fam wall $active{$qs}->[1]\n";
            $wall = substr($active{$qs}->[2],0,2);
            $strt = $active{$qs}->[0];
#            push @actfam, $fam if ($fam>0);
            $actfam{$fam} = $qs;  # a reverse hash
            $$rqs{$qs} = [$fam,$strt,$wall]  if ($fam>0);
    }
    return ($fndme,%actfam);
}       
sub testqs { #test if qsid still running 
    my $rqs = shift;
    my %active = ();
    my @stuff = ();
    my $qs;
    my $myjob = $ENV{"PBS_JOBID"};
    my $fndme = 0;
    $fndme = 1 if (! $myjob);   #non pbs system
    
#    print "testqs called @$rqs\n";
	$| = 1;    # flush STDOUT
#    print "enter testqs as $MYUSER\n";
    my $pid = open(QS,"qstat -e -u $MYUSER |");
    my $jnk = <QS>;
    $jnk = <QS>;
#    print "testqs $jnk\n";
    while (<QS>){
        chomp;
        @stuff = split;
        next if ($stuff[0] =~ /\D/);
#        print "qs $stuff[0], $stuff[-2]\n";
        if ($stuff[-2] =~ /R|Q/){
            $qs = "$stuff[0]";
            $fndme = 1 if ($qs == $myjob);
            $active{"$qs"} = 1;  #put found pids in hash
#            print "$stuff[0] is active\n";
        }
    }
    if (! $fndme) {
        print "did not find me job $myjob - something wrong - try later\n";
        return (0,0);
    }
    for  $qs (keys %$rqs){
#        print "test active $qs\n";
        next if ($qs =~ /D/);
#        print "$qs is still active\n" if exists $active{"$qs"};
        next if exists $active{"$qs"};
        my @blinfo = @{$$rqs{$qs}};

        print "$qs  blast $blinfo[0] not active\n" if ! exists $active{$qs};
        delete $$rqs{$qs};
        return ($qs,@blinfo);
    }
    return (0,0); #everthing still running or qued
}
