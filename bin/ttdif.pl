#!/usr/bin/env perl
use strict;
use Time::Piece;
use Time::Seconds;
use Cwd;
use Fcntl;
my $tw = time();
my $tstrt = localtime;
print "$tstrt now\n";
#open FW,">starttime";
#print FW "$tw\n";
#close FW;
sleep(3);
$tstrt = localtime;
print "$tstrt now\n";
if (-s 'starttime'){
    open FI,"<starttime";
    my $x = <FI>;
    close FI;
    chomp($x) if ($x);
    my $tw = time();
#    my $t = Time::Piece->strptime($x,"%+");
    my $dif = $tw - $x;
    my $tx = localtime($x);
    print "$tx in starttime, dif $dif \n";
}


