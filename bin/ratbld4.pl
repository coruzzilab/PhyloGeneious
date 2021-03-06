#!/usr/bin/env perl
## this is a paraphrase ofdameon's script to do a lot of rathchets
use Time::Piece;
use Time::Seconds;
use Cwd;
use Getopt::Long;
use POSIX qw/ceil/;
use strict;
my $t = localtime;
my %param = ( 'c' => 0,'h' =>1, 'i' => 1, 'm' => 'oid.nex', 'n' => 0, 'w' =>10, 'p' => 'pi', 't' => 24);
#getopts('gwp',\%param);
GetOptions("char:i" => \$param{'c'},
    "hold:i" => \$param{'h'},
    "iter:i" => \$param{'i'},
    "matrix:s" => \$param{'m'},
    "ncpu:i" => \$param{'n'},
    "prefix:s" =>\$param{'p'},
    "timeout:s" => \$param{'t'},
    "weight:i" =>\$param{'w'});

my $args = @ARGV; #number of args
my $c = $param{'c'};
my $h = $param{'h'};
my $i = $param{'i'};
my $m = $param{'m'};
my $n = $param{'n'};
my $w = $param{'w'};
my $prefix = $param{'p'};
print "h $h i $i m $m n $n \n";
my $timeout = $param{'t'};
my $timout = "";
my $timhr = $timeout;
my $timin = 0;
if ($timeout =~ /(\d+):(\d+):\d+/) {
    $timout = $timeout;
    $timhr = $1;
    $timin = $2;
}else{
    $timhr = $timeout;
    $timin = 0;
    $timout = "$timhr:00:00";
}
print "prefix will be $prefix timeout $timout\n";

if ($n and -s "$m") {   #mandatory parameters - should be in data subdirectory to run
    if ($h<9){
        $h = 9;
        print "hold must be 9 or more for pfijo\n";
    }

#	$w = int(($w/100)*$c);
	############################### MAKE PROC
	for(my $k = $n+1; $k > 1; $k--){
        my $prfil = "$prefix$k";
		my $buffer = "log +$prfil.log\n";
		$buffer .= "mxram 6000;\n";
		$buffer .= "taxname=;\n";
		$buffer .= "taxname +64;\n";
		$buffer .= "nstates 32;\n";
		$buffer .= "watch=;\n";
		$buffer .= "p $m;\n";
		$buffer .= 'hold ' . ($h + 1) . ";\n";
		$buffer .= 'rs ' . (int(rand(999999999))+1) . ";\n";
        $buffer .= "macro=;\n";
        $buffer .= "var: chk;\n";
        $buffer .= "set chk isfile $prefix$k.tre;\n";
        $buffer .= "if ('chk'==1) quote $prefix$k.tre exists exiting;zz;end;\n";
        $buffer .= "set chk isfile $prefix$k.rat;\n";
        $buffer .= "if ('chk'==1) short $prefix$k.rat;end;\n";
        $buffer .= "macro-;\n";

		$buffer .= "col 3;\n";
		$buffer .= "pfijo: slack 50 numsels 1;\n";
#		$buffer .= "cc [ .;\n";
# 		$buffer .= "cc /1 .;\n";
		$buffer .= "xi;\n";
		$buffer .= "riddup;\n";
		$buffer .= "mu=rep1ho" . $h . "wag;\n";
		$buffer .= "riddup-;\n";
		$buffer .= "tsave $prefix$k.rat;\n";
		$buffer .= "save.;\n";
		$buffer .= "timeout $timout;\n";
		$buffer .= "hold " . $h . ";\n";
		for(my $j = $i; $j > 0; $j--){############################### RAT
#			$buffer .= "quote ------------------ RATCHET WEIGHTS $j of $i;\n";
#			my @rat = ();
#            my @ctbl = ();
#            for (my $l=0;$l<=32;$l++){ #init table
#                $ctbl[$l] = ''; #start empty
#            }
#			for(my $x = $c; $x > 0; $x--){
#				$rat[$x] = 1;
#				}
#			for(my $x = $w; $x >= 0; $x--){
#				$rat[rand($c)] += 1;
#				}
#            for (my $x = $c;$x> 0;$x--){
#                my $v = $rat[$x];
#                $ctbl[$v] .= "$x " if $v<33;
#            }
#            for (my $l=1;$l<=32;$l++){
#                $buffer .= "cc /$l $ctbl[$l];\n" if $ctbl[$l];
#            }

			$buffer .= "quote ------------------ pinion (pfijo) SEARCH $j of $i;\n";
		    $buffer .= "xi;\n";
			$buffer .= "riddup;\n";
			$buffer .= "pfijo 0;\n";
			$buffer .= "riddup-;\n";
			$buffer .= "save.;\n";
			}
		$buffer .= "tsave/;\n";
		$buffer .= "k0;\n";
		$buffer .= "hold " . (($i * 3 * $h) + $h) . ";\n";
#		$buffer .= "cc [ .; cc /1 .; xi;\n";
		$buffer .= "short $prefix$k.rat;\n";
		$buffer .= "best;\n";
		$buffer .= "tsave $prefix$k.tre;\n";
		$buffer .= "save.;\n";
		$buffer .= "tsave/;\n";		
		$buffer .= "log/;\n";
		$buffer .= "zz;\n";
		$buffer .= "proc/;\n";
		open(OUTFILE, ">$prfil.proc") || die("could not open $prfil.proc\n!");
		print OUTFILE "$buffer";
		close OUTFILE
		}
        $h = $h*$n;
        my $prfil = "$prefix".'1';
        $timhr = $timhr*$n;
        $timin = $timin*$n;
        if ( $timin > 0) {
            $timhr += int($timin/60);
            $timhr = 8 if $timhr >8;
            $timin = $timin % 60;
            if ($timin>9) {
                $timout = "$timhr:$timin:00";
            }else{
                $timout = "$timhr:0$timin:00";
            }
        }else {
            $timhr = 8 if $timhr >8;
            $timout = "$timhr:00:00";
        }
		my $buffer = "log +$prfil.log\n";
		$buffer .= "mxram 6000;\n";
		$buffer .= "taxname=;\n";
		$buffer .= "taxname +64;\n";
		$buffer .= "nstates 32;\n";
		$buffer .= "watch=;\n";
		$buffer .= "report +180/10/10;\n";
		$buffer .= "p $m;\n";
		$buffer .= 'hold ' . ($h + 1) . ";\n";
		$buffer .= 'rs ' . (int(rand(999999999))+1) . ";\n";
        $buffer .= "macro [200;\n";
        $buffer .= "macro=;\n";
        $buffer .= "var: chk;\n";
    for (my $k = 2;$k<=$n+1;$k++){
        $buffer .= "set chk isfile $prefix$k.rat;\n";
        $buffer .= "if ('chk'==1) short $prefix$k.rat;end;\n";
    }
        $buffer .= "set chk ntrees+10;\n";
		$buffer .= "hold  'chk';\n";
        $buffer .= "macro-;\n";
        $buffer .= "tf: rounds 1 min 5;\n";
        $buffer .= "tsave $prfil.rat;\n";
		$buffer .= "timeout $timout;\n";
        $buffer .= "sort;\n";
        $buffer .= "score;\n";
		$buffer .= "xi;\n";
		$buffer .= "riddup;\n";
        $buffer .= "tf;\n";
		$buffer .= "riddup-;\n";
        $buffer .= "save .;\n";
		$buffer .= "riddup;\n";
        $buffer .= "tf;\n";
		$buffer .= "riddup-;\n";
        $buffer .= "save .;\n";
		$buffer .= "riddup;\n";
        $buffer .= "tf;\n";
		$buffer .= "riddup-;\n";
        $buffer .= "save .;\n";
		$buffer .= "riddup;\n";
        $buffer .= "tf;\n";
		$buffer .= "riddup-;\n";
        $buffer .= "save .;\n";
		$buffer .= "riddup;\n";
        $buffer .= "tf;\n";
		$buffer .= "riddup-;\n";
        $buffer .= "save .;\n";
        $buffer .= "tsave/ ;\n";
        $buffer .= "sort;\n";
        $buffer .= "score;\n";
        $buffer .= "best;\n";
        $buffer .= "score;\n";
        $buffer .= "tsave oids.tre;\n";
        $buffer .= "save .;\n";
        $buffer .= "tsave/ ;\n";
        $buffer .= "taxname ] ;\n";
        $buffer .= "tchoose 0 ;\n";
        $buffer .= "export - oid.tre;\n";
        $buffer .= "log/ ;\n";
        $buffer .= "system >.tre.done;\n";
		$buffer .= "zz;\n";
		$buffer .= "proc/;\n";
		open(OUTFILE, ">$prfil.proc") || die("could not open $prfil.proc\n!");
		print OUTFILE "$buffer";
		close OUTFILE


        
		
		

	} else { ############################### PRINT INSTRUCTIONS
		print("A Perl script for making TNT compatible pfijo procedure files\n");
		print("USAGE: ratbld4.pl  [-h hold] [-i iterations] -m matrix.ss -n files\n\n");
        print "[-t timeout (hours)]\n";
		print("WHERE\n");
#		print("\t-c\tnumber of characters in the matrix\n");
		print("\t-h\tnumber of trees to hold (default = $h)\n"); 
		print("\t-i\tnumber of ratchet iterations (default = $i)\n");
		print("\t-m\tmatrix file name\n");
		print("\t-n\tnumber of procedure files\n");
#		print("\t-w\tpercent to weight (default = $w)\n\n");
		print("\t-p\tfilename prefix (default = $prefix)\n\n");
		print("\t-t\ttimeout for a search (default = $timeout)\n\n");
		}



exit;


