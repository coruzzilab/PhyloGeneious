#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;
use Fcntl;
my $OID_HOME;
my $OID_USER_DIR;
my $OID_BLAST;
my $OID_BLASTDIR;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
    $OID_BLAST = "$OID_USER_DIR/blast";
}

use lib "$OID_HOME/lib";
use OrthologID;
use Getopt::Long;
use POSIX qw/ceil/;
use strict;
my $t     = localtime;
my @stuff = ();
my $strt;
my %qsid   = ();
my $outsfx = 1;
print "local time $t\n";

if ( @ARGV > 0 ) {
    $outsfx = "$ARGV[0]";
}
my $dummypart = 1000 + $outsfx;
my $reportfl  = "$OID_USER_DIR/mk_part$outsfx";

$OID_BLASTDIR = $OrthologID::OID_BLASTDIR;
my $comb  = $OrthologID::OID_COMB;
my $MAXQS = $OrthologID::MAXQS;
my $MINBL = $OrthologID::BLST_MIN;    #avg time in min of blast
my $secbl = $MINBL * 60;
print "I will try to run blasts at $secbl secs ($MINBL min) avg time\n";
print "blastdb $OID_BLASTDIR combined $comb tasks $MAXQS qsubs\n";
open( FA, "$comb" ) or die "could not find $comb\n";
my $nofa   = 0;
my $curfa  = "";
my @falens = ();
my $curlen = 0;
my $totlen = 0;
my $skip   = 100 * ( $outsfx - 1 );
print "we will skip $skip records\n";

while (<FA>) {
    chomp;
    if (/^>/) {
        if ($curfa) {
            push @falens, $curlen;
            $nofa++;
            $totlen += $curlen;
        }
        $curlen = length($_);
        $curfa  = substr( $_, 1 );    #skip >
        next;
    }
    $curlen += length($_);
}    #while
if ($curfa) {
    push @falens, $curlen;
    $nofa++;
    $totlen += $curlen;
}
die "no fa records\n" if ( $nofa < 0 );
## finised reading

print "there are $nofa genes, tot len $totlen\n";
my $fav   = $totlen / $nofa;
my $sumsq = 0;
foreach my $l (@falens) {
    $sumsq += ( $fav - $l ) * ( $fav - $l );
}
my $stdv = sqrt( $sumsq * 1. / $nofa );
print "avg $fav stdv $stdv\n";
die "not enough data to bother further\n " if ( $nofa < 100 );
## now get cute
close FA;
open( FB, ">$OID_BLAST/Part$dummypart.faa" )
  or die "couldn't open Part$dummypart.faa\n";
open( FA, "$comb" ) or die "could not find $comb\n";
my $rotlen = 0;
$curfa = '';
my $rused = 0;
my $use   = 100 * $fav;

while (<FA>) {
    chomp;
    if (/^>/) {
        if ( $skip > 0 ) {
            $skip--;
            next;
        }
        if ($curfa) {
            last if ( $rotlen > $use );
        }
        $curfa = substr( $_, 1 );    #skip >
        $rused++;
        print FB "$_\n";
        next;
    }
    $rotlen += length($_);
    print FB "$_\n";
}    #end while
print "output $rused fa records of len $rotlen\n";
close FA;
close FB;
$t = localtime;
print "submitting test blast at $t\n";
my $rc = system("rm $OID_BLAST/.?*$dummypart.done");
allBlast($dummypart);
my $ct = localtime;
print "blast done at $ct\n";
my $dif = $ct - $t;
$dif = 1 if ( $dif <= 0 );
my $forhr    = $secbl / $dif;
my $recforhr = ceil( $forhr * $rused );
my $sizforhr = $forhr * $rotlen;
print
"time $dif sec - need about $recforhr records/ len $sizforhr to run $MINBL min\n";
my $parts   = ceil( $totlen / $sizforhr );
my $siztous = $totlen / $parts;              #sizforhr leaves small last rec
my $dbs     = @OrthologID::INGROUP + @OrthologID::OUTGROUP;
my $thrash =
  ( $parts / $MAXQS ) * $MINBL;    # This is aprox total time with above calc
my $thrhr = $thrash / 60;

print
  "have $dbs differend blast dbs and need $parts  parts for $MINBL min runs\n";
print "compute total time about $thrash min or $thrhr hrs\n";
## if $secbl/$dbs < 300 may run too long
##       if ($dbs>12 and $parts < $dbs){
if ( $secbl / $dbs < 300 and $thrash > 60 ) {
    $forhr    = 300 * $dbs / $dif;
    $sizforhr = $forhr * $rotlen;
    $recforhr = ceil( $forhr * $rused );
    $parts    = ceil( $totlen / $sizforhr );

    print "above value would thrash cache as each blast would run ",
      $secbl / $dbs, " seconds\n";
    print "better: $parts parts: each runs for ", 5 * $dbs,
      " min : size $sizforhr, recs $recforhr\n";
    $thrash = ( $parts / $MAXQS * 300 * $dbs / 3600 );
    print "this vsn may do better than $thrash hrs\n";
}

#       $sizforhr = ceil($sizforhr);
$sizforhr = $siztous;
$rc       = system("rm $OID_BLAST/Part$dummypart.faa");
open( FB, ">$reportfl" ) or die "couldn't open $reportfl\n";
print FB "part\t$parts\tsize\t$sizforhr\trec\t$recforhr\tdif\t$dif\n";
close FB;
