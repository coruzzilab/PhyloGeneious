#!/usr/bin/env perl
my @sizes;
my $totsz = 0;
my $p = $ARGV[0];
$p = 1 if ! $p;
@sizes = `ls -s *Part$p.blst`;
foreach my $str(@sizes){
    my @st = split /\s/,$str;
    $totsz += $st[0];
}
print "size $totsz for part $p entries $#sizes\n";

