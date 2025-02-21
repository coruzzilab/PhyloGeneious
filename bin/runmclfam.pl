#!/usr/bin/env perl

use strict;
use Cwd;
use Time::Piece;

use Fcntl;
my $OID_HOME;
my $OID_USER_DIR;
my $OID_BLASTDB;
my $OID_BLASTDIR;
my $OID_DATADIR;
my $HPC;
my $ENV_WRAPPER;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $HPC          = $ENV{'HPC'};
    $HPC          = 'P' if !defined($HPC);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
    $OID_BLASTDIR = "$OID_USER_DIR/blast";
    $OID_BLASTDB  = "$OID_USER_DIR/blastdb";
    $OID_DATADIR  = "$OID_USER_DIR/data";
    $ENV_WRAPPER  = $ENV{'ENV_WRAPPER'};
    if (!defined($ENV_WRAPPER)){
        $ENV_WRAPPER = "";
    }
}

use lib "$OID_HOME/lib";
use OrthologID;
my @pidtbl  = ();
my @opentbl = ();
my @clustbl = ();
our @lastfam = ();
my %fastast = ();
my $nclust;
my $nthr = 2;
my $tm;

my $mcl_exe      = "mcl";
my $clustersFile = "$OID_BLASTDIR/clusters";
my $mciFile      = "$OID_BLASTDIR/weights.mci";
my $tabFile      = "$OID_BLASTDIR/weights.tab";
my $clusterdone  = "$OID_BLASTDIR/.mci.done";
my $mcxdone      = "$OID_BLASTDIR/.mcx.done";
my $familydone   = "$OID_DATADIR/.family.done";
my $singletFile  = "$OID_DATADIR/singlets";
$tm = localtime;
my $savDir = getcwd;
print "runmclfam.pl in $savDir at $tm\n";

#    chdir "1";
#    die "runtntmx.pl nproc nthreads prefix\n" if (@ARGV<3);
die "runmclfam.pl ncpu\n" if ( @ARGV < 1 );
my $nthr = $ARGV[0];
print "we run on $nthr cpus\n";
if ( -f $familydone ) {
    print "family already made - nothing to do\n";
    exit(0);
}
if ( !-f $clusterdone ) {    #must run mci
                             # run mcl
    print "$tm Running mcl ...\n";
    my @mclArgs = (
        $mciFile,   "-force-connected", "y",  "-scheme", "7", "-I", "1.4",
        "-use-tab", $tabFile,           "-o", $clustersFile
    );
    push( @mclArgs, "-te", $nthr );
    if (! $ENV_WRAPPER eq "") {
        my $status = system($ENV_WRAPPER, $mcl_exe, @mclArgs);
        die "$mcl_exe exit status = $?" unless $status == 0;
    }
    else {
        my $status = system( $mcl_exe, @mclArgs );
        die "$mcl_exe exit status = $?" unless $status == 0;
    }
    $tm = localtime;
    print "$tm mcl done.\n";
    open( XX, ">$clusterdone" );    #create a done file
    close XX;
}
## now finished with mcl - read in clusters file
open MFH, $clustersFile or die "Cannot open $clustersFile: $!\n";
while (<MFH>) {
    chomp;
    push @clustbl, $_;
}

# we read in entire set of genes into memory - why use blastcmd
my $genect = 0;
my $filect = 0;
for my $faa (<$OID_BLASTDB/*.faa>) {

    # i don't trust combined.fa -
    my $geneid;
    $filect++;
    open FA, $faa or print "cannot open $faa ";
    while (<FA>) {
        if (/^>(\S+)/) {
            $geneid = $1;
            $fastast{$1} = "";    #new gene
            $genect++;
        }
        else {
            $fastast{$geneid} .= $_;    #append whole record
        }
    }    # read each faa file
    close FA;
}    # all faa files
print "we loaded $genect genes from $filect files\n";

$nclust = scalar(@clustbl);
for ( my $i = 0 ; $i < $nthr ; $i++ ) {
    $lastfam[$i] = 0;
}
## now we can start to make families
$tm = localtime;
print "$tm starting to make $nclust families\n";
for ( my $th ; $th < $nthr ; $th++ ) {
    my $pid;

    #        sleep 2;
    #        next if $pid = fork;  #parent
    next if ( $pid = fork() );    # parent

    # now in child
    for ( my $fam = $th ; $fam < $nclust ; $fam += $nthr )
    {                             #we chose families based on thread
        my $sizfm = make1fam( $fam + 1, $clustbl[$fam], \%fastast );
        if ( $sizfm < 2 ) {
            last;
        }
    }

    #child just completed
    exit(0);
}    #end of threads

1 while ( wait() != -1 );
$tm = localtime;
print "$tm all threads done\n";

$savDir = getcwd;
chdir $OID_DATADIR;

my $count = 0;
my %seen;
foreach my $dir (<[1-9]*>) {
    $seen{$dir} = 1;
    $count++;
}
chdir $savDir;
print "found $count directories\n";
my @genes = split /\s/, $clustbl[$count];
if ( @genes != 1 ) {    # implies something missing

    for ( my $i = 1 ; $i <= $nclust ; $i++ ) {
        next if ( $seen{$i} );
        @genes = split /\s/, $clustbl[ $i - 1 ];
        last if ( @genes == 1 );

        print "missing dir $i try again\n";
        $count++;
        chdir $savDir;
        my $sizfm = make1fam( $i, $clustbl[ $i - 1 ], \%fastast );
        print "make1fam returns $sizfm\n";
    }
}    #if missing directory

# now do singlets
open SFH, ">>$singletFile"
  or die "Create open $singletFile for writing: $!\n";
my $singletNum = 0;
for ( my $i = $count ; $i < $nclust ; $i++ ) {
    @genes = split /\s/, $clustbl[$i];
    if ( @genes == 1 ) {
        print SFH "$genes[0]\n";
        $singletNum++;
    }
    else {
        my $j = $i + 1;
        print "bug found fam $j where expected singlet genes (@genes)\n";
        my $sizefm = make1fam( $i + 1, $clustbl[$i], \%fastast );
    }
}
$tm = localtime;
print "$tm made $count families $singletNum singlets\n";

open( XX, ">$familydone" );    #create a done file
close XX;
