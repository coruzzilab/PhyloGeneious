#!/usr/bin/env perl
package Qsblast . pm;
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

our @EXPORT = qw(runBlast
  )

  my $OID_CONF = "blastcfg";
my $Blastall = "blastp";
my $blstdir  = "";
our $num_thr = 16;
our $MAXQS   = 10;
our $mem_bl  = "16GB";
our $eval    = "1e-10";
our @dbfiles = ();
our $COMBIN  = "comb.fa";
our $PARM    = "";
our $COMBOUT = "comb.blst";
my $infile;
my $cmdline = "";
our $myjob;
our $HPC = "S";

# Initialize
&initblstid;

sub init_myjob {

    if ( $HPC =~ /S/ ) {
        $myjob = $ENV{"SLURM_JOB_ID"};
        print "my job $myjob\n" if ($myjob);
    }
    else {
        $myjob = $ENV{"PBS_JOBID"};
        print "my job $myjob\n" if ($myjob);
    }
}

sub initblstid {
    open CONF_FH, "$OID_CONF" or die "Cannot open $OID_CONF: $!\n";
    while (<CONF_FH>) {
        chomp;
        next if (/^#/);
        if (/^\s*HPC\s*=\s*([PSC]).*$/) {
            $HPC = $1;
        }
        elsif (/^\s*BLAST\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
            $Blastall = $1;
        }

        elsif (/^\s*DB\s*=\s*(.+?)\s*$/) {       # Ingroup species prefixes
            $_       = $1;
            @dbfiles = split;
        }
        elsif (/^\s*COMBIN\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
            $COMBIN = $1;
        }
        elsif (/^\s*COMBOUT\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
            $COMBOUT = $1;
        }
        elsif (/^\s*PARM\s*=\s*(.+?)\s*$/) {       #any other parameter
            $PARM = $1;
        }
        elsif (/^\s*MEM\s*=\s*(\d+)/) {            # Number of CPUs available
            $mem_bl = "$1GB";
        }
        elsif (/^\s*EVAL\s*=\s*(.+?)\s*$/) {       # Ingroup species prefixes
            $eval = $1;
        }
        elsif (/^\s*NTHR\s*=\s*(\d+)/) {           # Number of CPUs available
            $num_thr = $1;
        }
        elsif (/^\s*VERBOSE\s*=\s*(\d+)/) {        # Verbosity
            $VERBOSE = $1;
        }
        elsif (/^\s*DEBUG\s*=\s*(\d+)/) {          # Debug mode
            $DEBUG = $1;
        }
        elsif (/^\s*BLSTMIN\s*=\s*(\d+)/) { #avg time in minute to set for blast
            $BLST_MIN = $1;
        }

        elsif (/^\s*MAXQS\s*=\s*(\d+)/) {    # Maximum qsub jobs
            $MAXQS = $1;
        }
    }
    close CONF_FH;
    init_myjob();
}    #end of init`

sub runBlast {
    my $prefix = shift;

    my $blstinfile = "$blastdir/Part$prefix.faa";    #instead use original
    my ( $blastCmd, @blastOut );

    print "BLASTing .. $prefix\n";
    foreach my $sp (@dbfiles) {
        my $t = localtime;
        print "$t db sb $sp\n";
        my $outfile = "$blastdir/$sp.Part$prefix.blst";
        if ( system( "grep -q '\\.$sp\\.Part$prefix\\.done' $blastdir/.$sp.Part.done") < 1 ) {
            print "$sp.Part$prefix done - skipping\n" if $verbose > 1;
            next;
        }

# Run blast
#		$blastCmd =
#		    "$BLASTALL -a $numProc -p blastp -d $sp -I -e $eCutOff -m 8 -i $tmpFile";
        $blastCmd =
            "$Blastall -db $sp -query $blstinfile -evalue $eval "
          . "-num_threads $numProc  -out $outfile $PARAM";

        my $blastOutput = `$blastCmd`;
        $t = localtime;
        print "$t complete blast $sp rc $blastOutput\n" if $verbose > 1;
        open( XX, ">$blastdir/.$sp.Part.done" );    #add to done file
            print XX ".$sp.Part$prefix.done\n";
        close XX;
    }
}    ## runBlast
