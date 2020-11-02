#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;
use Fcntl;
use AnyDBM_File;
use DB_File;
my $OID_HOME;
my $OID_USER_DIR;
my $OID_DATADIR;
my $NCPU;   # cpu's to give blast
our $MAXQS; 
our $TNTB;  #threshold for running tnt job (less we skip)
my $gb;   #gb to run blast
my $MYUSER;
my $HPC;
my $myjob;
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);
    $OID_DATADIR="$OID_USER_DIR/data";
    $MYUSER = `whoami`;
}

use lib "$OID_HOME/lib";
use OrthologID; 
use Getopt::Long;
use strict;
my %sztaxa = ();
my %ncfam = ();
my %rcdir = (); # attempt to recover
my $Checkpoint="$OID_USER_DIR/log/job/chkbigmon.log";
my $Chkfa;
my $Chkcnt=0;
my %chktbl = ();
my %maftbl = ();

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

sub xfmthpc{
    (my $gb, my $wall,  my $cpus, my $tasks, my $script, my @args) = @_;
    if  ($HPC =~ /S/){
        my $PARAMS = "-n $tasks -c $cpus -t $wall --mem $gb ";
        my $argstr = '';
        foreach my $arg(@args){
            $argstr .= " $arg";
        }
        my $submit = q/sbatch /."$PARAMS $script $argstr";
        return $submit;
    } else {   #HPC = 'P' or undefined
       my $walmem = "mem=$gb".",walltime="."$wall";
       my $PARAMS = "-l nodes=1:ppn=$NCPU -l $walmem";
       my $argstr = "";
       if (@args >0){
            $argstr = "-v";
            my $argnum = 0;
            foreach my $arg (@args){
              $argnum++;
              $argstr .= " arg$argnum"." $arg"; 
          }
        }   
        my $submit =  q/qsub  /."$PARAMS $script $argstr";
        return ($submit);
    }
}
sub xgetqids {
    my $rqs = shift;
    my $flg = 1;
    my @stuff = ();
    my $qs;

    if  ($HPC =~ /S/){
        my $pid = open(QS,"squeue -l -u $MYUSER |");
        my $jnk = <QS>;
        $jnk = <QS>;
        while (<QS>){
            chomp;
            @stuff = split;
            if ($flg) {
#                print " squeue " ,join " ",@stuff;
                print "squeue ",scalar(@stuff);
#                print "\nX$stuff[0]X\n";
#                $flg = 0;
            }
            next if ($stuff[0] =~ /\D/);
            $qs = "$stuff[0]";
            $$rqs{"$qs"} = 1; #mark active
            print "$qs and $$rqs{$qs}\n";
        }

        return;
    } else {

        my $pid = open(QS,"qstat -e -u $MYUSER |");
        my $jnk = <QS>;
        $jnk = <QS>;
#    print "testqs $jnk\n";
        while (<QS>){
            chomp;
            @stuff = split;
            next if ($stuff[0] =~ /\D/);
            if ($stuff[-2] =~ /R|Q/){
                $qs = "$stuff[0]";
                $$rqs{"$qs"} = 1; #mark active
            }
        }
    }
}

sub startqs { #starts a qsub
#                startqs(\@qsid,$cmpqs,$famx);
    my ($rqs,$fam,$noftsk,$thr,$prefix,$wall);
    my ($JOB_SCRIPT,$maft);
    if (@_>3){
        ($rqs, ,$fam, $noftsk, $thr, $prefix, $wall ) = @_;
        $JOB_SCRIPT="$OID_HOME/bin/run_tntmx".$HPC.".sh";
        $maft = 0;
    } else {
        ($rqs, $fam) = @_;
        print "startqs for mafft for $fam\n";
        $JOB_SCRIPT="$OID_HOME/bin/run_mafft".$HPC.".sh";
        $maft = 1;
    }
        
#    my $walmem = OrthologID::setlimit($sztaxa{$fam});
    my ($PARAMS,$mem,$submit);
    my $nt = 1;
    my %gb;
    if (! $maft) {
       $mem = 4*$thr;
       $mem += $mem if $sztaxa{$fam} > 5000;
       $mem = 64 if $noftsk == 4; #very big 
#    my $walmem = "mem=$mem"."GB,walltime=$wall"; 
#    my $JOB_SCRIPT="$OID_USER_DIR/run_tntmx.sh";
#       $PARAMS = "-l nodes=1:ppn=$thr,walltime=$wall -l mem=$mem"."GB ";
#      $submit = q/qsub  /." $PARAMS $JOB_SCRIPT -v arg1=$fam,arg2=$noftsk,arg3=$thr,arg4=$prefix 2>/dev/null";
      $gb = "$mem"."GB";
      my @args = ();
      push @args,$fam;
      push @args,$noftsk;
      push @args,$thr;
      push @args,$prefix;

      $submit = fmthpc($gb,$wall,$thr,$nt,$JOB_SCRIPT,@args);
      $submit .= " 2>/dev/null";
    }else {
        if ($sztaxa{$fam}<4000){
            $mem = "16GB";
            $thr = 2;
            $noftsk = -1;
            $wall = "04:00:00";
        }else{
            $mem = "32GB";
            $thr=4;
            $noftsk = -1;
            $wall = "08:00:00";
        }
#      $PARAMS = "-l nodes=1:ppn=2,walltime=04:00:00 -l mem=$mem ";
#      $submit = q/qsub /." $PARAMS $JOB_SCRIPT -v arg1=$fam,arg2=2 2>/dev/null";
#      $wall = "04:00:00";
      $prefix="mafft";
      my @args = ();
      push @args,$fam;
      push @args,$thr;
      $submit = fmthpc($mem,$wall,$thr,$nt,$JOB_SCRIPT,@args);
      $submit .= " 2>/dev/null";
  }
    print "$submit\n";
#    my $pid = open(QS,"qsub -l $PARAMS $JOB_SCRIPT -v arg1=/"$at/",arg2=/"^$fam.$dol/" 2>/dev/null |");
    my $pid = open(QS,"$submit |");
    
    my $qs = <QS>;
    if ($qs){
        chomp $qs;
        $qs = fixqs($qs);
        print "start qsid $qs for family $fam\n";
        $$rqs{$qs} = $fam; #put $qs in its place
        my $tm = localtime;
        if (exists $rcdir{$fam}) {
            my $retry = $rcdir{$fam}->[0]+1;
            $rcdir{$fam} = [$retry,$tm,$wall];
        }else{
            $rcdir{$fam} = [0,$tm,$wall];
        }
        if($maft){
            $maftbl{$fam} = 1;
        }else{
            delete $maftbl{$fam} if exists $maftbl{$fam};
        }
        $chktbl{$qs} = [$fam,$noftsk,$thr,$prefix,$wall];
        $Chkcnt++;
        if ($Chkcnt<50){
            my $line = "$qs\t$fam\t$noftsk\t$thr\t$prefix\t$wall\n";
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
    my $fam;
    my %actfam = ();
#    my $myjob = $ENV{"PBS_JOBID"};
    my $fndme = 0;
    my $flag = 1;
    $fndme = 1 if (! $myjob);   #non pbs system
#    print "testqs called @$rqs\n";
	$| = 1;    # flush STDOUT
#    print "enter testqs as $MYUSER\n";
    getqids(\%active);
    if ($myjob){
        $fndme = 1 if exists $active{"$myjob"};
    }
    if (! $fndme) {
        print "did not find me job $myjob - something wrong - try later\n";
        return ($fndme,%actfam);
    }
    foreach $qs (keys %active){
        if (exists $chktbl{$qs}){
            ($fam,@stuff) = @{$chktbl{$qs}};
            print "active $qs found fam $fam wall $active{$qs}\n";
            if ($fam>0){
                if ($stuff[0] < 0) { #this is a maft job
                    $maftbl{$fam} = 1;
                }
                $actfam{$fam} = $qs;  #reverse hash
                $$rqs{$qs} = $fam;
            }
            next;
        }

        
#        my $pid = open(GR, "grep -c '>' FAMILY |") ;
#        $sztaxa{$dir} = <GR>;
        my $log = "log/job/*".$qs."*";
#        print $log
        my $pid = open(GR, "grep 'runtntmx ' $log 2>/dev/null| ");
        $fam = 0;
        $_ = <GR>;
        next if (! $_);
        chomp;
#        print "found $_\n";
        if (/data\/(\d+)/){
            $fam = $1;
        }
        print "active $qs found fam $fam wall $active{$qs} by log\n";
        $actfam{$fam} = $qs if ($fam>0);  #reverse hash
        $$rqs{$qs} = $fam  if ($fam>0);
    }
    return ($fndme,%actfam);
}       
sub setupoid{
    my $setup_proc = "$OID_HOME/bin/setup_proc.sh";
    my $id = shift;
    my $ntax = $sztaxa{$id};
    my $nchar = $ncfam{$id};
    print "$id has $nchar  char and $ntax taxa \n";
    my $rc = system("$setup_proc $id $ntax $nchar");
    open (FA,"setup.txt");
    $_ = <FA>;
    chomp;
    my @stuff = split;
    close FA;
    return(@stuff);

}
sub chkfam{
    my $rqs = shift;
    my $fam = shift;
    my $oldqs = shift;
    my $famdir = "$OID_DATADIR/$fam";
    if (-f "$famdir/.tre.done"){
        delete $rcdir{$fam};
        return 1; # completed ok
    }
    my $tm = localtime;

    if (not exists $maftbl{$fam}) {
        if (! exists $rcdir{$fam}) {
            print "chkfam $fam no record\n";
            $rcdir{$fam} = [0,$tm,"24:00:00"];
        }

        my @oldv = @{$rcdir{$fam}};
        print "$fam $tm old $oldv[1] wall $oldv[2]\n";
        if ($oldv[0]>1) {
            print "$fam retry failed $oldv[0] times\n";
            return 1; #give up
        }

# new code - .nnnnnn.start has starttime
        $tm = time();
        my $sttm;
        if (open(QS,".$oldqs.start") ) {
            $sttm = <QS>;
            close QS;
        }else {
            print ".$oldqs.start not found \n";
            $sttm = $tm;
        }
        my $dif = $tm-$sttm;
        my $wallv = substr $oldv[2],0,2;
        $wallv = $wallv*3600;
        my $savdir = getcwd;
        chdir "$OID_DATADIR/$fam";
        my @stuff = setupoid($fam);
        chdir $savdir;
        if ($dif < $wallv-600) { #then did not time out
            print "old only ran $dif sec - no timeout just retry\n";
            my $rc =  startqs($rqs,$fam,@stuff);
            return 0 if $rc;
            print "startup fails for $fam\n";
            return 1;
        }
        $wallv = substr $oldv[2],0,2;
        $wallv = int($wallv*3/2);
        my $oldwall = pop  @stuff;
        substr($oldwall,0,2) = $wallv; #up time
#    $oldwall = substr $oldwall,0,2,$wallv; #up time
        push @stuff, $oldwall;
        print "increase time to $oldwall\n";
        my $rc =  startqs($rqs,$fam,@stuff);
        return 0 if $rc;
        print "startup fails for $fam\n";
        return 1;
    }else {   # completed mafft
        delete $maftbl{$fam};
        delete $rcdir{$fam} if exists $rcdir{$fam};
        my $savdir = getcwd;
        chdir "$OID_DATADIR/$fam";
        if (! -s "oid.nex") {
            print "mafft failed for $fam\n";
            chdir $savdir;
            return 1;
        }
        print "completed mafft for $fam\n";
            open (NX, "oid.nex");
            $_=<NX>;
            $_=<NX>;
            close NX;
            chomp;
            my @stuff = split;
            if ($stuff[1] != $sztaxa{$fam}) {
                print "oid.nex for $fam has nchar $stuff[0], ntaxa $stuff[1]\n";
                print "not correspond to my db\n";
            }
            $ncfam{$fam} = $stuff[0];
            @stuff = setupoid($fam);
            chdir $savdir;
            my $rc =  startqs($rqs,$fam,@stuff);
            return 0 if $rc;
            print "startup fails for $fam\n";
            return 1;
        }
        
}


### main program
my $t = localtime;
my $strtm = localtime; ## for restart
my $qsubed;
my %qsid = ();
my $qsct = 0;
my $oldqs;
my $qsfam = 0;
my $schedrst = 0;
my %actfam = ();
my $procnt = 0;
our $familyDB = "$OID_USER_DIR/familyDB.dbm";
    $MAXQS = $OrthologID::MAXQS;
    $TNTB = $OrthologID::TNTB;
    print "maxqs $MAXQS\n";
print "local time $t\n";
##our ($opt_g, $opt_w, $opt_p);
    my $savdir = getcwd;
    my $badtie = 0;
    tie (%sztaxa,"DB_File",$familyDB,O_RDONLY,0666) or
        do {print "could not create tie $familyDB\n";$badtie = 1;};
    if ($badtie) {
 	 chdir $OID_DATADIR;
     foreach my $dir (<[1-9]*>) {
         chdir $dir;
         my $pid = open(GR, "grep -c '>' FAMILY |") ;
         $_ = <GR>;
         chomp;
         $sztaxa{$dir} = $_;
 	    chdir $OID_DATADIR;
      }
     }
     chdir $OID_USER_DIR;
    my @kefnd = keys %sztaxa;
    my $fam = scalar @kefnd;
    print "there are $fam families\n";
#    my @taxsz = sort { $sztaxa{$a} <=> $sztaxa{$b} } @kefnd;
    my @taxsz = sort { $a <=> $b} @kefnd;
    $HPC = $OrthologID::HPC; #env nogood
    print "we are in hpc $HPC\n";
    if ($HPC =~ /S/){
        $myjob = $ENV{"SLURM_JOB_ID"};
        print "my job $myjob\n";
    }else {
        $myjob = $ENV{"PBS_JOBID"};
        print "my job $myjob\n";
    }
        rdchkpt();  #reads in checkpoint file if it exists, else creates it

    my $fndme = 0;
    until ($fndme) {  ## sometimes qstat doesn't work - so wont find me
        ($fndme,%actfam) = activqs(\%qsid);
        sleep 15 if (! $fndme);
    }
    writechkpt(\%qsid);
        $qsct = keys %qsid;  #init with actual active tasks
        print "were running $qsct tasks on restart\n";
            while ($qsct >= $MAXQS ){
                ($oldqs,$qsfam) = OrthologID::testqs(\%qsid);
                delete $chktbl{$oldqs} if ($oldqs and exists $chktbl{$oldqs});
                if ($oldqs and chkfam(\%qsid,$qsfam,$oldqs)) {
                    print "completed family $qsfam qcnt $qsct\n";
                    $qsct-- ;
                    last;
                }
                sleep(60);
            }  #end while
  for (my $redo=0;$redo <2;$redo++){
    for $fam (@taxsz) {
#        print "looking at $fam\n";
        chdir "$OID_DATADIR/$fam";
        next if (-f ".tre.done"); #finished
        next if ($sztaxa{$fam} < $TNTB);  #could use last, but this is safer if strange order
        next if ($actfam{$fam});
        next if ($rcdir{$fam}); #already running  - if we do a second pass important
        my @parms = ();
        if (-s "oid.nex"){
            open (NX, "oid.nex");
            $_=<NX>;
            $_=<NX>;
            close NX;
            chomp;
            my @stuff = split;
            if ($stuff[1] != $sztaxa{$fam}) {
                print "oid.nex for $fam has nchar $stuff[0], ntaxa $stuff[1]\n";
                print "not correspond to my db\n";
                next;
            }
            $ncfam{$fam} = $stuff[0];
            @parms = setupoid($fam);
        } else { #must run mafft for fam
            @parms = ();
        }
            
            chdir $OID_USER_DIR;
            $qsct++ if (startqs(\%qsid,$fam,@parms));
            while ($qsct >= $MAXQS ){
                print "waiting $qsct after family $fam\n";
                ($oldqs,$qsfam) = OrthologID::testqs(\%qsid);
                delete $chktbl{$oldqs} if ($oldqs and exists $chktbl{$oldqs});
                if ($oldqs and chkfam(\%qsid,$qsfam,$oldqs)) {
                    print "completed family $qsfam qcnt $qsct\n";
                    delete $actfam{$qsfam} if exists $actfam{$qsfam};
                    $qsct-- ;
                    last;
                }
                sleep(60);
            }  #end while
    }
  }
    chdir $savdir;
    my $slept=0;
            while ($qsct > 0 ){
                ($oldqs,$qsfam) = OrthologID::testqs(\%qsid);
                delete $chktbl{$oldqs} if ($oldqs and exists $chktbl{$oldqs});
                if ($oldqs and chkfam(\%qsid,$qsfam,$oldqs)) {
                    print "completed family $qsfam qcnt $qsct\n";
                    $qsct-- ;
                    next;
                }
                next if ($oldqs);
                sleep(60);
                $slept++;
                if ($slept>5){
                    print "wainting for jobs to complete\n";
                    $slept=0;
                }
            }  #end while
            open (QS,">log/job/bigmdone");  #tell schedular
            print "completed all families\n";

