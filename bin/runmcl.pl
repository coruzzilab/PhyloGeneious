#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;
use Fcntl;
my $OID_HOME;
my $OID_USER_DIR;
my $OID_BLASTDIR;
my $OID_DATADIR;
our $NCPU;   # cpu's to give blast
our $glsttm;  # external shell script start time in sec
my $ncpub = 0;
my $MYUSER;
my $OID_MCL;
my $myjob;
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);
    $OID_BLASTDIR="$OID_USER_DIR/blast";
    $OID_DATADIR="$OID_USER_DIR/data";
    $MYUSER = `whoami`;
    $OID_MCL=$ENV{'OID_MCL'};
    $OID_MCL=2 if ! defined($OID_MCL);
    print "using $OID_MCL GB per thread for mcl\n";
}

use lib "$OID_HOME/lib";
use OrthologID; 
use Getopt::Long;
use strict;
our $HPC;
my $clusterdone = "$OID_BLASTDIR/.mci.done";
my $mcxdone = "$OID_BLASTDIR/.mcx.done";
my $mcifile = "$OID_BLASTDIR/weights.mci";
my $familydone = "$OID_DATADIR/.family.done";
my $BLAST_RES_DB   = "$OID_BLASTDIR/blastres.blst";

my $Checkpoint="$OID_USER_DIR/log/job/chkmclfam.log";
my $Chkfa;
my $Chkcnt=0;
my %chktbl = ();
my $sleeptm = 300;

sub rdchkpt{   #read in checkpoint file if it exists
    %chktbl = ();
    if (-s "$Checkpoint"){
        open (FA,"$Checkpoint") or do {print "failed open of checkpoint $Checkpoint\n";
            open (FA,">$Checkpoint"); }
    }else{
        open (FA,">$Checkpoint");
    }
    while (<FA>){
        chomp;
        next if (!$_); 
        my @chk = split;
        my $pid = @chk[0];
        $chktbl{$pid} = [@chk[1 .. $#chk]];
    }
    close FA;
}
sub writechkpt{  #rewrite entire chkpoint
    my $rqs = shift;   #%qsid has actual active
    my $cnt = 0;
    my $fa;
    open (FA,">$Checkpoint");
    foreach my $pid (keys %chktbl){
        if (! exists $$rqs{$pid}) {  #only use active pids
            delete $chktbl{$pid};
            next;
        }
        $cnt++;
        my $line = "$pid\t";
        my @chk = @{$chktbl{$pid}};
        foreach my $item (@chk){
            $line .= "\t$item";
        }
        print FA "$line\n";
    }
    print "output $cnt lines to $Checkpoint\n";
    close FA;
    open ($fa,">>$Checkpoint");
    $Chkfa = $fa;
    $Chkcnt=0;
}




sub startqs { #starts a qsub
#                startqs(\@qsid,$blid,$wall);
    (my $JOB_SCRIPT, my $rqs, my $gbx, my $blid, my $ncpux, my $wl, my $parg) = @_;
    #    my $JOB_SCRIPT="$OID_USER_DIR/run_oid_job.sh";

#    my $submit = q/qsub  /."$PARAMS $JOB_SCRIPT".q/ -v arg1="-b"/.",arg2=$blid  2>/dev/null";
    my $nt = 1;
    my $wall = "$wl:00:00";
    my $gb = "$gbx"."GB";
    my $submit = fmthpc($gb,$wall,$ncpux,$nt,$JOB_SCRIPT,@$parg);
    $submit .= " 2>/dev/null";
    print "$submit\n";
#    my $pid = open(QS,"qsub -l $PARAMS $JOB_SCRIPT -v arg1=/"$at/",arg2=/"^$fam.$dol/" 2>/dev/null |");
    my $pid = open(QS,"$submit |");
    my $t = localtime;
    my $qs = <QS>;
    if ($qs){
        chomp $qs;
        $qs = fixqs($qs);
        print "start qsid $qs for blastid $blid\n";
#        if ($qs !~ /^\d+$/){
#            my @stf = split /\s/,$qs;
#            if ($stf[-1] =~ /^\d+$/){
#                $qs = $stf[-1];
#            } else {
#                print "truble breaking up $qs\n";
#            }
#        }

        my $wallhr = substr($wall,0,2);
#        $$rqs{$qs} = [$blid,$wallhr]; #put $qs in its place
        $$rqs{$qs} = $blid; #put $qs in its place
        $chktbl{$qs} = [$blid,$wallhr];
        $Chkcnt++;
        if ($Chkcnt<50){
            my $line = "$qs\t$blid\t$wallhr\n";
            print $Chkfa "$line"; 
        }else{
            writechkpt($rqs);
            $Chkcnt=0;
        }
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
    my $blid;
    my $wallhr;
    my %actfam = ();
    my $strt;
    my $wall;
    my $fndme = 0;
    my $flag = 1;
    $fndme = 1 if (! $myjob);   #non pbs system
    my $t = localtime;
#    print "testqs called @$rqs\n";
	$| = 1;    # flush STDOUT
    getqids(\%active);
    if ($myjob){
        $fndme = 1 if exists $active{"$myjob"};
    }
#        return ($fndme,%actfam) if ($fndme);

#    print "enter testqs as $MYUSER\n";
    foreach $qs(keys %active){
        next if($qs == $myjob);
        if (not exists $chktbl{"$qs"}) {
            print "$qs active but not in checkpoint\n";
            delete $active{"$qs"};
            next;
        }
        ($blid,$wallhr) = @ {$chktbl{"$qs"}};
#        $$rqs{"$qs"} = $chktbl{"$qs"};
        $$rqs{"$qs"} = $blid;
        $actfam{"$blid"} = $qs;
        print "found $qs running for blast $blid\n";
    }
    foreach $qs(keys %chktbl){
        delete $chktbl{"$qs"} if (not exists $active{"$qs"});
    }
    return ($fndme,%actfam);
}       

sub oldtestqs { #test if qsid still running 
    my $rqs = shift;
    my %active = ();
    my @stuff = ();
    my $qs;
    my $fndme = 0;
    $fndme = 1 if (! $myjob);   #non pbs system
    
#    print "testqs called @$rqs\n";
	$| = 1;    # flush STDOUT
    getqids(\%active);
#    my @kq = keys %active;
#    print "keys active ",join ",",@kq,"\n";
#    print "enter testqs as $MYUSER\n";
    if ($myjob){
        $fndme = 1 if exists $active{"$myjob"};
    }
    if (! $fndme) {
        print "did not find me job $myjob - something wrong - try later\n";
        return (0,0);
    }
    for  $qs (keys %$rqs){
#        print "test active $qs\n";
        next if ($qs =~ /D/);
#        print "$qs is still active\n" if exists $active{"$qs"};
        next if exists $active{$qs};
#        foreach my $x (%active){
#            print "$x not number\n " if ($x =~/D/);
#            print "$x and $active{$x}\n";
#            print "$qs is a keys\n " if ($qs == $x);
#        }
        print "really keys active ",join ":",keys %active,"\n";
#        my @blinfo = @{$$rqs{$qs}};

        my $blinfo = $$rqs{$qs};
#        print "$qs  blast $blinfo not active\n" if ! exists $active{$qs};
        delete $$rqs{$qs};
#        return ($qs,@blinfo);
        return ($qs,$blinfo);
    }
    return (0,0); #everthing still running or qued
}

sub waitqs{
    my $pqsid = shift;
    while(1){
        my $oldqs;
        my $qsfam;
        ($oldqs,$qsfam) = OrthologID::testqs($pqsid);
        delete $chktbl{$oldqs} if ($oldqs and exists $chktbl{$oldqs});
        return if($oldqs);
        sleep $sleeptm;
        chkuptm();
    }
}
### main program
my $t = localtime;
my $qsubed;
if (-f $familydone) {
    print "$t family done - nothing to do\n";
    exit(0);
}
$HPC = $OrthologID::HPC; #env nogood
print "$t we are in hpc $HPC\n";
if ($HPC =~ /S/){
    $myjob = $ENV{"SLURM_JOB_ID"};
    print "my job $myjob\n";
}else {
    $myjob = $ENV{"PBS_JOBID"};
    print "my job $myjob\n";
}
my $mcx_script = "$OID_HOME/bin/run_mcx$HPC.sh";
my $mclfam_script = "$OID_HOME/bin/run_mclfam$HPC.sh";
my $nthr = 0;
my $memavail=0;

my $args = @ARGV; #number of args
die "runmcl availmem (ncpu)\n" if($args < 1); 
my $memarg = $ARGV[0];
$memarg =~ /(\d+)GB/;
$memavail = $1;
die "need at least 12GB = got $memarg ($memavail)\n" if $memavail < 12;
$nthr = $ARGV[1] if ($args>1);
$nthr = $OrthologID::NCPU if ($nthr <1);
print "runmcl will assume $memavail GB and $nthr thread for mcl\n";

   

my $qsct = 0;
my %qsid = ();
my %actblst = ();
my $tskid;
my $bldone;

### startup code
rdchkpt(); # this is last run (unless done)
$qsubed = 0;
until($qsubed){
        ($qsubed,%actblst) = activqs(\%qsid);
        sleep 15 if (! $qsubed);
    }
    writechkpt(\%qsid);
         
    $qsct = keys %qsid;  #init with actual active tasks
    my @kq = keys %qsid;
    print "were running $qsct tasks on restart ",join ":",@kq,"\n";
    setglsttm(); #get checkpointed time
    my $strtm = $OrthologID::glsttm;
    if ($qsct) {
        $bldone = waitqs(\%qsid);
    }
        die "cannot find blastres.blst\n" if (! -f  $BLAST_RES_DB);
        my @blstar = `ls -s $BLAST_RES_DB`;
        my $blstsz = $blstar[0];
        $blstsz =~ s/^\s//;
        my @blstrec = split /\s/,$blstsz;
        my $gb = int($blstrec[0]/1000000);
        $gb = 1 if !$gb;
        print "blast size $gb GB($blstrec[0])\n";
        my $needmem = int($gb/3); #approx mem for mcxload is 1/3 of blastres size
        my $needtm = (($gb+24)/25) * 3600 ; # needs about 1hr per 25GB
        $needmem = 1 if ($needmem <1);
        $needtm = 1 if ($needtm < 1);
        $sleeptm = 30 if $gb < 3; # small blast quick run
    while ( ! -f $mcxdone ) { # need to do mcxload
        my  $tmnw = time();
        my $td = $tmnw - $strtm;
        if($memavail >= $needmem and (($td + $needtm) < 72000)){
            print "need $needtm time and $needmem GB have $memavail gb and enough time so run mcx inline\n";
            mclmcx();
        }else {
            print "need to scedule mcxload\n";
            my @args;
            my $wl = int(($gb+24)/25);
            $wl = 8 if($wl < 8);
            $qsubed = 0;
            until($qsubed){
                if(startqs($mcx_script,\%qsid,$needmem,1,1,$wl,\@args)) {
                    $qsubed = 1;
                }
            }
            $bldone = waitqs(\%qsid);
        }
    }
    while ( ! -f $familydone ) {
        my $needmem = $nthr * $OID_MCL;
        my $wl = 16 ; #rediculus - shoud be  < 2 hrs
        my @args;
        push @args,$nthr;
        $qsubed = 0;
        until($qsubed){
            if (startqs($mclfam_script,\%qsid,$needmem,2,$nthr,$wl,\@args)){
               $qsubed = 1;
            }
        }
        $bldone = waitqs(\%qsid);
    }
    my $tm = localtime;
    print "$tm runmcl done\n";    
