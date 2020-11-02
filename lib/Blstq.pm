#####!/usr/bin/env perl
package Blstq;
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
	'all' => [
		qw(

		  )
	]
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
    runBlast
    fixqs
    fmthpc
    getqids
    combBlast
);

my $OID_HOME;        # OrthologID main directory
my $OID_USER_DIR;    # OrthologID user run directory
die "OID_HOME environment variable undefined!\n" if !$ENV{'OID_HOME'};
die "OID_USER_DIR environment variable undefined!\n" if !$ENV{'OID_USER_DIR'};
$OID_HOME       = $ENV{'OID_HOME'};
$OID_USER_DIR   = $ENV{'OID_USER_DIR'};
my $OID_CONF="blastcfg";
my $Blastall="blastp";
our $Blastdir="";
my $Dbdir = "";
our $num_thr=16;
our $NCPU;            # Number of CPU to be used
our $MAXQS=10;
our $mem_bl="16GB";
our $eval="1e-10";
our @dbfiles=();
our $COMBIN="comb.fa";
our $PARAM="";
our $COMBOUT="comb.blst";
our $BLST_MIN=60;
my $VERBOSE=0;
my $DEBUG=0;
my $infile;
my $cmdline="";
our $myjob;
our $HPC="S";
my $MYUSER = `whoami`;
# Initialize
&initblstid;

sub init_myjob{
    
    if ($HPC =~ /S/){
        $myjob = $ENV{"SLURM_JOB_ID"};
        print "my job $myjob\n" if ($myjob);
    }else {
        $myjob = $ENV{"PBS_JOBID"};
        print "my job $myjob\n" if ($myjob);
    }
}

sub initblstid{
	open CONF_FH, "$OID_CONF" or die "Cannot open $OID_CONF: $!\n";
	while (<CONF_FH>) {
		chomp;
		next if (/^#/);
        if (/^\s*HPC\s*=\s*([PSC]).*$/) 
        {
            $HPC = $1;
        }
		elsif (/^\s*BLAST\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$Blastall       = $1;
        }
		elsif (/^\s*BLASTDIR\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$Blastdir       = $1;
        }
		elsif (/^\s*DBDIR\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$Dbdir       = $1;
        }

		elsif (/^\s*DB\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$_       = $1;
			@dbfiles = split;
		}
		elsif (/^\s*COMBIN\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$COMBIN       = $1;
        }
		elsif (/^\s*COMBOUT\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$COMBOUT       = $1;
        }
        elsif (/^\s*PARAM\s*=\s*(.+?)\s*$/) {   #any other parameter
            $PARAM = $1;
        }
		elsif (/^\s*MEM\s*=\s*(\d+)/) {            # Number of CPUs available
			$mem_bl = "$1GB";
		}
		elsif (/^\s*EVAL\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$eval = $1;
		}
		elsif (/^\s*NTHR\s*=\s*(\d+)/) {            # Number of CPUs available
			$num_thr = $1;
            $NCPU=$num_thr;
		}
		elsif (/^\s*VERBOSE\s*=\s*(\d+)/) {         # Verbosity
			$VERBOSE = $1;
		}
		elsif (/^\s*DEBUG\s*=\s*(\d+)/) {           # Debug mode
			$DEBUG = $1;
		}
        elsif (/^\s*BLSTMIN\s*=\s*(\d+)/) {           #avg time in minute to set for blast
            $BLST_MIN = $1;
        }

		elsif (/^\s*MAXQS\s*=\s*(\d+)/) {           # Maximum qsub jobs
			$MAXQS = $1;
		}
	}
	close CONF_FH;
    init_myjob();
    print "Blstq for $Blastall in dir $Blastdir\n";
    return 1;
}  #end of init`

sub fixqs {
    my $qs = shift;
        if ($qs !~ /^\d+$/){
            my @stf = split /\s/,$qs;
            if ($stf[-1] =~ /^\d+$/){
                $qs = $stf[-1];
                return $qs;
            } else {
                print "truble breaking up $qs\n";
                return 0;
            }
        }
        return $qs;
    } #end fixqs

sub fmthpc{
    (my $gb, my $wall,  my $cpus, my $tasks, my $script, my @args) = @_;
    if  ($HPC =~ /S/){
        my $PARAMS = "-N 1 --ntasks-per-node $tasks -c $cpus -t $wall --mem $gb ";
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
sub runBlast {
    my $prefix = shift;

    my $blstinfile = "$Blastdir/Part$prefix.faa";  #instead use original
	my ( $blastCmd, @blastOut );
    my $verbose = $VERBOSE;
    my $t;
    

	print "BLASTing .. $prefix\n" ;
	foreach my $sp ( @dbfiles ) {
        my $sppth = "$Dbdir/$sp";
        $t = localtime;
        print "$t db sb $sp\n";
        my $outfile = "$Blastdir/$sp.Part$prefix.blst";
        if ( -f "$Blastdir/.$sp.Part$prefix.done"){
            print "$sp.Part$prefix done - skipping\n" if $verbose>1; 
            next;
        }

		# Run blast
#		$blastCmd =
#		    "$BLASTALL -a $numProc -p blastp -d $sp -I -e $eCutOff -m 8 -i $tmpFile";
            $blastCmd = 
                "$Blastall -db $sppth -query $blstinfile -evalue $eval "
                . "-num_threads $num_thr  -out $outfile $PARAM";
            print "command: $blastCmd\n";
                 
		my $blastOutput = `$blastCmd`;
        $t = localtime;
        print "$t complete blast $sp rc $blastOutput\n" if $verbose>1;
        open(XX,">$Blastdir/.$sp.Part$prefix.done");  #create a done file
        close XX;
    }
        print "$t completed all blast for Part $prefix\n" if $verbose>1;
        open(XX,">$Blastdir/.$prefix.done");  #create a done file
        close XX;
}  ## runBlast

sub combBlast {
    my $noparts = shift;
    my $cntrec = 0;
    open (FB, ">$Blastdir/$COMBOUT") or do { print "open for $COMBOUT fails\n";return };
    for (my $part=1;$part <= $noparts;$part++) {
        foreach my $sp (@dbfiles) {
            open FA,"$Blastdir/$sp.Part$part.blst" or do {print "no $sp.Part$part.blast found.. skipping \n";next};
            while (<FA>){
                $cntrec++;
                print FB $_;
            }
            close FA;
        } #foreach $sp
    } #for parts
    close FB;
    print "copied $cntrec records to $Blastdir/$COMBOUT\n";
}
