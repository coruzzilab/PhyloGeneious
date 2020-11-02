#!/usr/bin/env perl
use strict;
my $locus;
my %fanm;
my $dupcnt = 0;
print $ARGV[0],"\n";
while (<>){
    if (/^>(\S+)/){
        $locus = $1;
        if (exists $fanm{$locus}){
            print "dupe locus ",$locus,"\n";
            $dupcnt++;
        } else{
            $fanm{$locus} = 1;
        }
    }
}
print $dupcnt, "dupes\n";
