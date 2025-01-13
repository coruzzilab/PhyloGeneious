#!/usr/bin/env perl
use strict;
use Cwd;
my $locus;
my %fanm;
my $dupcnt  = 0;
my %ntaxa   = ();
my @thisdis = ( 0, 0, 0 );
print $ARGV[0], "\n";
my $savDir = getcwd;

foreach my $dir (<[1-9]*>) {
    chdir $dir;
    my $pid = open( GR, "grep -c '>' FAMILY |" );
    $ntaxa{$dir} = <GR>;
    chdir $savDir;
}
my @kefnd = keys %ntaxa;
my $fam   = scalar @kefnd;
print "there are $fam families\n";
my @taxsz = sort { $ntaxa{$a} <=> $ntaxa{$b} } @kefnd;
for $fam (@taxsz) {
    my $thissz = $ntaxa{$fam};
    if ( $thissz < 100 ) {
        $thisdis[0]++;
    }
    elsif ( $thissz < 500 ) {
        $thisdis[1]++;
    }
    else {
        $thisdis[2]++;
        print "$fam has size $thissz\n";
    }
}
print "counts @thisdis\n";

