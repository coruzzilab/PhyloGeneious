#!/usr/bin/env perl
my @sizes;
my $totsz = 0;
my $avsz = 0;
my $maxsz = 0;
if ( -f "blastres.blst") {
    my @blstar = `ls -s blastres.blst`;
    my $blstsz = $blstar[0];
    $blstsz =~ s/^\s//;
    my @blstrec = split /\s/,$blstsz;
    my $mb = $blstrec[0]/1000;
    print "blast size $mb ($blstrec[0])\n";
}
if (! -f "Part1.faa") {
    print "Part1 missing\n";
    exit(1);
}
my $Parts = `ls Part*faa | wc -l`;
for (my $p=1;$p<=$Parts;$p++){
    @sizes = `ls -s *Part$p.blst`;
    $totsz = 0;
    foreach my $str(@sizes){
        $str =~ s/^\s//;
         my @st = split /\s/,$str;
        $totsz += $st[0];
    }
    print "size $totsz for part $p \n";
    $avsz += $totsz;
    $maxsz = $totsz if $totsz>$maxsz;
}
$avsz = $avsz/$Parts;
print "average of $Parts is $avsz and biggest $maxsz\n";

