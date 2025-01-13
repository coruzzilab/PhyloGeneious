#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;
use Fcntl;
my $OID_HOME;
my $OID_USER_DIR;
my $Blastdir;
my $HPC;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $HPC          = $ENV{'HPC'};
    $HPC          = 'P' if !defined($HPC);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
}

use lib "$OID_HOME/lib";
use Blstq;
use Getopt::Long;
use POSIX qw/ceil/;
use strict;
my $t     = localtime;
my @stuff = ();
my $strt;
my $outsfx = 1;
print "local time $t\n";
my $dummypart = 1000 + $outsfx;
my $parts     = 0;
my $ptsz      = 0;
my $ptrec     = 0;
my $ptdif     = 0;
my $bigused   = 0;
my $runs      = 2;
my $reportfl  = "$OID_USER_DIR/mk_partq$outsfx";

$Blastdir = $Blstq::Blastdir;
my $comb  = $Blstq::COMBIN;
my $MAXQS = $Blstq::MAXQS;
my $MINBL = $Blstq::BLST_MIN;    #avg time in min of blast
my $TCPU  = $Blstq::num_thr;     #cpu's to start for blast
my $MEM   = $Blstq::mem_bl;      #amount of memory to give
$HPC = $Blstq::HPC;              #env nogood
print "$comb $MAXQS $MINBL $TCPU from blastq\n";
my $secbl = $MINBL * 60;
print "I will try to run blasts at $secbl secs ($MINBL min) avg time\n";

if ( -f "$Blastdir/.Part.done" ) {
    print "$Blastdir/.Part.done exists...\n";
    print "parts files already made -- skipping\n";
    exit(0);
}
print "blastdb $Blastdir combined $comb tasks $MAXQS qsubs\n";
## first run 2 dummy blasts to compute sizes
my $run_shell = "$OID_HOME/bin/run_mk_blastq" . $HPC . ".sh";
print "$run_shell\n";
for ( my $i = 1 ; $i <= $runs ; $i++ ) {
    my $reportflx = "$OID_USER_DIR/mk_partq$i";
    unlink $reportflx;

    #           system ("$OID_USER_DIR/qstest.sh $i");
    #           system ("$OID_HOME/bin/run_mk_blast.sh $i $TCPU $MEM");
    my $rc = system("$run_shell $i $TCPU $MEM");
    print "rc from system $rc\n";

    #           if ($rc != 0){
    #               die "system bad rc $rc\n";
    #           }
}
## now spin and wait for the jobs to complete
my $fndfile = $runs;
print "looking for $reportfl\n";
sleep 20;
while ($fndfile) {
    if ( -s $reportfl ) {
        $fndfile--;
        if ( $fndfile > 0 ) {
            $outsfx++;
            $reportfl = "$OID_USER_DIR/mk_partq$outsfx";
            print "looking for $reportfl\n";
            next;    #all could be completed by now
        }
        else {
            last;    #just to save time
        }

    }
    sleep 20;
}
print "found all parts\n";
for ( $outsfx = 1 ; $outsfx <= $runs ; $outsfx++ ) {
    $reportfl = "$OID_USER_DIR/mk_partq$outsfx";
    open( FA, "$reportfl" ) or die "cannot open $reportfl\n";
    $_ = <FA>;
    chomp;
    @stuff = split /\s/;

# now stuff[1] has parts, stuff [3] has approx size of part stuff[5] has approx no fa's in part
# stuff[7] is timediff for this run
    print "$reportfl: ", join "\t", @stuff, "\n";
    $parts += $stuff[1];
    $ptsz  += $stuff[3];
    $ptrec += $stuff[5];
    $ptdif += $stuff[7];
}

# calulate average
$parts = ceil( $parts / $runs );
$ptsz  = ceil( $ptsz / $runs );
$ptrec = ceil( $ptrec / $runs );
$ptdif = $ptdif / $runs;
print "we will use part $parts size $ptsz (recs ~ $ptrec, time $ptdif)\n";
open( FA, "$comb" ) or die "could not find $comb\n";
my $curfa = '';    #needed between parts
my @fatxt = ();

for ( my $prefix = 1 ; $prefix <= $parts ; $prefix++ ) {
    open( FB, ">$Blastdir/Part$prefix.faa" )
      or die "couldn't open Part$dummypart.faa\n";
    my $rotlen = 0;
    my $rused  = 0;

    while (<FA>) {
        chomp;
        if (/^>/) {
            if ($curfa) {
                if ( !@fatxt
                    or ( $fatxt[0] =~ /^X+$/ and $fatxt[$#fatxt] =~ /^X+$/ ) )
                {
                    print "bad fa entry $curfa\n";
                }
                else {
                    print FB ">$curfa\n";
                    foreach my $fa (@fatxt) {
                        print FB "$fa\n";
                        $rotlen += length($fa);
                    }
                }
            }
            @fatxt = ();                 #must empty
            $curfa = substr( $_, 1 );    #skip >
            $rused++;
            last if ( $rotlen >= $ptsz );
            next;
        }
        push @fatxt, $_;

        #           $rotlen += length($_);
    }    #end while
    if ( !$_ and $curfa ) {    #ot eof
        if ( !@fatxt or $fatxt[0] =~ /^X+$/ ) {
            print "bad fa entry $curfa\n";
        }
        else {
            print FB ">$curfa\n";
            foreach my $fa (@fatxt) {
                print FB "$fa\n";
                $rotlen += length($fa);
            }
        }
    }
    print "PART$prefix.faa output $rused fa records of len $rotlen\n";
    $bigused += $rotlen;
    close FB;
}
close FA;
$bigused = $bigused / $parts;
print "made $parts partfiles ave len $bigused\n";
open( FA, ">$Blastdir/.Part.done" )
  ;    # mark sucessful completion of making parts so it isn't redone
close FA;
