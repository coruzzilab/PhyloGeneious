#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;

#sleep(60);
my $t    = localtime;
my $dayv = 3600 * 24;    #day in sec
for (@ARGV) {
    my $tv = -A;
    $tv = $tv * $dayv;    #convert to seconds
    my $fst = $t - $tv;
    print "$_ created $fst\n";
}
