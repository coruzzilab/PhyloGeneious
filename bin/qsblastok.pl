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
my $myjob;
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);
    $OID_BLAST="$OID_USER_DIR/blast";
    $MYUSER = `whoami`;
}

use lib "$OID_HOME/lib";
use OrthologID; 
use Getopt::Long;
use strict;
our $HPC;

my $Checkpoint="$OID_USER_DIR/log/job/chkqsblast.log";
my $Chkfa;
my $Chkcnt=0;
my %chktbl = ();

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
sub fmthpc{
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
       my $narg = @args;
       if (@args >0){
            $argstr = "-v ";
            my $argnum = 0;
            foreach my $arg (@args){
              $argnum++;
              $argstr .= "arg$argnum"."=$arg"; 
              $argstr .= "," if ($argnum<$narg);
          }
        }   
        my $submit =  q/qsub  /."$PARAMS $script $argstr";
        return ($submit);
    }
}
sub getqids {
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
#                startqs(\@qsid,$blid,$wall);
    (my $rqs, ,my $blid, my $wall) = @_;
    my $JOB_SCRIPT="$OID_USER_DIR/run_oid_job.sh";

#    my $submit = q/qsub  /."$PARAMS $JOB_SCRIPT".q/ -v arg1="-b"/.",arg2=$blid  2>/dev/null";
    my $nt = 1;
    my @args = ();
    push @args,'"-b"';
    push @args,$blid;
    my $submit = fmthpc($gb,$wall,$NCPU,$nt,$JOB_SCRIPT,@args);
    $submit .= " 2>/dev/null";
    print "$submit\n";
#    my $pid = open(QS,"qsub -l $PARAMS $JOB_SCRIPT -v arg1=/"$at/",arg2=/"^$fam.$dol/" 2>/dev/null |");
    my $pid = open(QS,"$submit |");
    my $t = localtime;
    my $qs = <QS>;
    if ($qs){
        chomp $qs;
        print "start qsid $qs for blastid $blid\n";
        if ($qs !~ /^\d+$/){
            my @stf = split /\s/,$qs;
            if ($stf[-1] =~ /^\d+$/){
                $qs = $stf[-1];
            } else {
                print "truble breaking up $qs\n";
            }
        }

        my $wallhr = substr($wall,0,2);
        $$rqs{$qs} = [$blid,$wallhr]; #put $qs in its place
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
        $$rqs{"$qs"} = $chktbl{"$qs"};
        $actfam{"$blid"} = $qs;
        print "found $qs running for blast $blid\n";
    }
    foreach $qs(keys %chktbl){
        delete $chktbl{"$qs"} if (not exists $active{"$qs"});
    }
    return ($fndme,%actfam);
}       
sub testqs { #test if qsid still running 
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
        my @blinfo = @{$$rqs{$qs}};

        print "$qs  blast $blinfo[0] not active\n" if ! exists $active{$qs};
        delete $$rqs{$qs};
        return ($qs,@blinfo);
    }
    return (0,0); #everthing still running or qued
}
### main program
my $t = localtime;
my $strtm = localtime; ## for restart
my $qsubed;
print "local time $t\n";
$HPC = $OrthologID::HPC; #env nogood
print "we are in hpc $HPC\n";
if ($HPC =~ /S/){
    $myjob = $ENV{"SLURM_JOB_ID"};
    print "my job $myjob\n";
}else {
    $myjob = $ENV{"PBS_JOBID"};
    print "my job $myjob\n";
}


##our ($opt_g, $opt_w, $opt_p);
my %param = ( 'g' => 16,'w' =>24, 'q' => 20, 'n' => 10);
my $help='';
GetOptions("g:i" => \$param{'g'},
    "help|?" => \$help,
    "w:i" => \$param{'w'},
    "n:i" => \$param{'n'},
    "q:s" => \$param{'q'});

die "qsblast.pl noblastfiles (-g nn -w nn -q nn -n cpus) gb,wall and max qsubs cpus for blast\n" if $help;
my $args = @ARGV; #number of args
my $fir = 0 ; #first positional arg
print "args $args , def  GB(g) $param{'g'}, wall(w) $param{'w'} maxqs(q) $param{'q'} ncpu(n) $param{'n'}\n";

for (my $i=0;$i<$args;$i++){
    print "$i,$ARGV[$i]\n";
        }
my $srch = 'Part.*faa';
my $nofbl =0;
if ($args>0) {  ##
    $nofbl = int($ARGV[0]);
}
    
my $srch = 'Part.*faa';
if ($nofbl<1) {
#         my $pid = open(GR, qw/grep -c "$srch" $OID_BLAST/." 2>/dev/null | ");
       my $grst = "ls $OID_BLAST | grep -c '$srch'  ";
       print "$grst\n";
       my $pid = open(GR,"$grst |");
         $nofbl = <GR>;
         chomp $nofbl;
     } ## scan for number of part files to blast

#my $nofbl = $ARGV[0];
print "will do $nofbl blasts\n";
$gb = "$param{'g'}"."GB";
my $wall = "$param{'w'}".":00:00";
my $MAXQS = "$param{'q'}";
$NCPU = $param{'n'};  
my $blastp = 0;
my $qsct = 0;
my %qsid = ();
my %actblst = ();
#my $lastblast = 0;

### startup code
rdchkpt(); # this is last run (unless done)
$qsubed = 0;
until($qsubed){
        ($qsubed,%actblst) = activqs(\%qsid);
        sleep 15 if (! $qsubed);
    }
    writechkpt(\%qsid);
         
#        if (@actblst){
#            @actblst = sort {$a <=> $b} @actblst;
#            $lastblast = $actblst[$#actblst];
#        }
        $qsct = keys %qsid;  #init with actual active tasks
        print "were running $qsct tasks on restart \n";
for (my $bl=1;$bl<= $nofbl;$bl++){
    $blastp++;
    print "look for $OID_BLAST/.$bl.done\n";
    next if ( -f "$OID_BLAST/.$bl.done" );  # did this guy aready
    next if ( $actblst{$bl});
    if ($qsct >= $MAXQS) {
        $blastp--;  #we never started up last part
        print "reaches MAXQS $MAXQS procs $qsct at blast $blastp\n";
        $|; #flush 
        last;
    }
    $qsubed = 0;
    until($qsubed){
        if (startqs(\%qsid,$bl,$wall)) {
            $qsct++;
            $qsubed = 1;
        }else {
            print "qsub fails at blast $blastp wait a while\n";
            sleep 60;
        }
    }
}
print "last blast submitted $blastp\n" if ($blastp >= $nofbl);

while (1){
    my $oldqs;
    my $bl;
    my $waltm;
    my $blastdone=0;
    my $subtm;
    last if ($qsct <=0); # nothing to wait for
    ($oldqs,$bl,$waltm) = testqs(\%qsid);
    if (! $oldqs) {
        sleep 60;
        my $tmnw = localtime;
        my $td = $tmnw - $strtm;  #when we start
        if($td > 72000) {
            print "time now $tmnw, diff $td\n";
            OrthologID::restartme();
            }
        next;
    }
    print "$oldqs,$bl,$waltm from testqs\n";
    delete $actblst{$bl} if (exists $actblst{$bl});
    $t = localtime;
    $qsct--;
    if ( -f "$OID_BLAST/.$bl.done"){
        print "blast $bl done at $t\n";
        next if ($blastdone);
       while (1) {  #needed because we may have to traverse busy blasts
        $blastp++;
        if ($blastp> $nofbl){
            print "last blast $nofbl was submitted\n";
            $blastdone = 1;
            last;
        }
        $bl = $blastp;
        next if ( -f "$OID_BLAST/.$bl.done");
        next if ( $actblst{$bl});
        $qsubed = 0;
        until($qsubed){
            if (startqs(\%qsid,$bl,$wall)) {
                $qsct++;
                $qsubed = 1;
            }else {
                print "qsub fails at blast $blastp wait a while\n";
                sleep 300;
            }
        }
        last; #we submitted job - exit while
      }  #internal while to find next blast
    }else { # something wront with blast
# new code - .nnnnnn.start has starttime
        $t = time();
        my $sttm;
        if (open(QS,".$oldqs.start") ) {
            $sttm = <QS>;
            close QS;
        }else {
            print ".$oldqs.start not found \n";
            $sttm = $t;
        }
        my $dif = $t - $sttm;
        my $walsec = $waltm * 3600 - 600; # give a leway of 10min
       if ($dif < $walsec) {
            print "blast canceled for $bl but clock not run out\n";
       } else {
        print "gave blast $bl wall $waltm but timed out\n";
        $waltm = int ($waltm*1.5); #convert back to int
        print "try new wall time $waltm\n";
       }  # don't die - maybe i crashed 
       $waltm= "$waltm:00:00";
        $qsubed = 0;
        until($qsubed){
            if (startqs(\%qsid,$bl,$waltm)) {
                $qsct++;
                $qsubed = 1;
            }else {
                print "qsub fails at blast $blastp wait a while\n";
                sleep 300;
            }
        }
    } #something wrong with blast
} # loop until all tasks done
$t = localtime;
print "$t all blasts done\n";
my  $rc = system("$OID_HOME/bin/orthologid.pl -B");
$t = localtime;
print "$t ran ortholodid -B rc $rc\n";
