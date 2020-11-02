#!/usr/bin/env perl
## this is a paraphrase ofdameon's script to do a lot of rathchets
use Time::Piece;
use Time::Seconds;
use Cwd;
use Getopt::Long;
use POSIX qw/ceil/;
use strict;
my $t = localtime;
my %param = ( 'c' => 0,'h' =>1, 'i' => 200, 'm' => 'oid.nex', 'n' => 0, 'w' =>10);
#getopts('gwp',\%param);
GetOptions("char:i" => \$param{'c'},
    "hold:i" => \$param{'h'},
    "iter:i" => \$param{'i'},
    "matrix:s" => \$param{'m'},
    "ncpu:i" => \$param{'n'},
    "weight:i" =>\$param{'w'});

my $args = @ARGV; #number of args
my $c = $param{'c'};
my $h = $param{'h'};
my $i = $param{'i'};
my $m = $param{'m'};
my $n = $param{'n'};
my $w = $param{'w'};
print "c $c h $h i $i m $m n $n w $w\n";

if ($c and $n and -s "$m") {   #mandatory parameters - should be in data subdirectory to run

	$w = int(($w/100)*$c);
	############################### MAKE PROC
	for(my $k = $n; $k > 0; $k--){
        my $prfil = $k+1;
		my $buffer = "log +$prfil.log\n";
		$buffer .= "mxram 4000;\n";
		$buffer .= "taxname=;\n";
		$buffer .= "taxname +64;\n";
		$buffer .= "nstates 32;\n";
		$buffer .= "p $m;\n";
		$buffer .= 'hold ' . ($h + 1) . ";\n";
		$buffer .= 'rs ' . (int(rand(999999999))+1) . ";\n";
        $buffer .= "macro=;\n";
        $buffer .= "var: chk;\n";
        $buffer .= "set chk isfile $k.rat;\n";
        $buffer .= "if ('chk'==1) short $k.rat;end;\n";
        $buffer .= "macro-;\n";

		$buffer .= "col 3;\n";
		$buffer .= "cc [ .;\n";
		$buffer .= "cc /1 .;\n";
		$buffer .= "xi;\n";
		$buffer .= "riddup;\n";
		$buffer .= "mu=rep1ho" . $h . "TBR;\n";
		$buffer .= "tsave $k.rat+;\n";
		$buffer .= "riddup-;\n";
		$buffer .= "save.;\n";
		$buffer .= "hold " . $h . ";\n";
		for(my $j = $i; $j > 0; $j--){############################### RAT
			$buffer .= "quote ------------------ RATCHET WEIGHTS $j of $i;\n";
			my @rat = ();
			for(my $x = $c; $x > 0; $x--){
				$rat[$x] = 1;
				}
			for(my $x = $w; $x >= 0; $x--){
				$rat[rand($c)] += 1;
				}
			my $one = ();
			my $two = ();
			my $three = ();
			my $four = ();
			my $five = ();
			my $six = ();
			my $seven = ();
			my $eight = ();
			my $nine = ();
			for(my $x = $c; $x >= 0; $x--){
				if($rat[$x] == 1){
					$one .= "$x ";
					} elsif($rat[$x] == 2){
						$two .= "$x ";
						} elsif($rat[$x] == 3){
							$three .= "$x ";
							} elsif($rat[$x] == 4){
								$four .= "$x ";
								} elsif($rat[$x] == 5){
									$five .= "$x ";
									} elsif($rat[$x] == 6){
										$six .= "$x ";
										} elsif($rat[$x] == 7){
											$seven .= "$x ";
											} elsif($rat[$x] == 8){
												$eight .= "$x ";
												} elsif($rat[$x] == 9){
													$nine .= "$j ";
													}	
				}
			if(length($one)){
				$buffer .= "cc /1 $one;\n";
				}
			if(length($two)){
				$buffer .= "cc /2 $two;\n";
				}
			if(length($three)){
				$buffer .= "cc /3 $three;\n";
				}
			if(length($four)){
				$buffer .= "cc /4 $four;\n";
				}	
			if(length($five)){
				$buffer .= "cc /5 $five;\n";
				}
			if(length($six)){
				$buffer .= "cc /6 $six;\n";
				}
			if(length($seven)){
				$buffer .= "cc /7 $seven;\n";
				}
			if(length($eight)){
				$buffer .= "cc /8 $eight;\n";
				}
			if(length($nine)){
				$buffer .= "cc /9 $nine;\n";
				}
			$buffer .= "quote ------------------ RATCHET SEARCH $j of $i;\n";
			$buffer .= "riddup;\n";
			$buffer .= "bbreak=tbr;\n";
			$buffer .= "riddup-;\n";
			$buffer .= "save.;\n";
			$buffer .= "quote ------------------ JACKKNIFE DEACTIVATIONS;\n";
			$buffer .= "cc /1 .;\n";
			$buffer .= "cc ]";
			for(my $x = $c; $x >= 0; $x--){
				if((int(rand(9999))+1) < 3679){
					$buffer .= " " . $x;
					}
				}
			$buffer .= ";\n";
			$buffer .= "riddup;\n";
			$buffer .= "quote ------------------ JACKKNIFE SEARCH $j of $i;\n";
			$buffer .= "bbreak=tbr;\n";
			$buffer .= "riddup-;\n";
			$buffer .= "save.;\n";	
			$buffer .= "quote ------------------ EQUAL WEIGHT SEARCH $j of $i;\n";
			$buffer .= "cc [ .;\n";
			$buffer .= "xi;\n";
			$buffer .= "riddup;\n";
			$buffer .= "bbreak=tbr;\n";
			$buffer .= "riddup-;\n";
			$buffer .= "save.;\n";
		    $buffer .= "tsave/ ;\n";
		    $buffer .= "tsave $k.rat+;\n";
			}
		$buffer .= "tsave/;\n";
		$buffer .= "k0;\n";
		$buffer .= "hold " . (($i * 3 * $h) + $h) . ";\n";
		$buffer .= "cc [ .; cc /1 .; xi;\n";
		$buffer .= "short $k.rat;\n";
		$buffer .= "best;\n";
		$buffer .= "tsave $k.tre;\n";
		$buffer .= "save.;\n";
		$buffer .= "tsave/;\n";		
		$buffer .= "log/;\n";
		$buffer .= "zz;\n";
		$buffer .= "proc/;\n";
		open(OUTFILE, ">oi$prfil.proc") || die("could not open oi$prfil.proc\n!");
		print OUTFILE "$buffer";
		close OUTFILE
		}
		my $buffer = "log +1.log\n";
		$buffer .= "mxram 5000;\n";
		$buffer .= "taxname=;\n";
		$buffer .= "taxname +64;\n";
		$buffer .= "nstates 32;\n";
		$buffer .= "p $m;\n";
		$buffer .= 'hold ' . ($h + 1) . ";\n";
		$buffer .= 'rs ' . (int(rand(999999999))+1) . ";\n";
        $buffer .= "macro [200;\n";
        $buffer .= "macro=;\n";
        $buffer .= "var: chk;\n";
    for (my $k = 1;$k<=$n;$k++){
        $buffer .= "set chk isfile $k.rat;\n";
        $buffer .= "if ('chk'==1) short $k.rat;end;\n";
    }
        $buffer .= "macro-;\n";
        $buffer .= "sort;\n";
        $buffer .= "score;\n";
        $buffer .= "tf;\n";
        $buffer .= "sort;\n";
        $buffer .= "score;\n";
        $buffer .= "best;\n";
        $buffer .= "score;\n";
        $buffer .= "tsave oi1.tre;\n";
        $buffer .= "save .;\n";
        $buffer .= "tsave/ ;\n";
        $buffer .= "log/ ;\n";
		$buffer .= "zz;\n";
		$buffer .= "proc/;\n";
		open(OUTFILE, ">oi1.proc") || die("could not open oi1.proc\n!");
		print OUTFILE "$buffer";
		close OUTFILE


        
		
		

	} else { ############################### PRINT INSTRUCTIONS
		print("A Perl script for making TNT compatible ratchet procedure files\n");
		print("USAGE: ratchet.pl -c characters [-h hold] [-i iterations] -m matrix.ss -n files\n[-w \%weight]\n\n");
		print("WHERE\n");
		print("\t-c\tnumber of characters in the matrix\n");
		print("\t-h\tnumber of trees to hold (default = $h)\n"); 
		print("\t-i\tnumber of ratchet iterations (default = $i)\n");
		print("\t-m\tmatrix file name\n");
		print("\t-n\tnumber of procedure files\n");
		print("\t-w\tpercent to weight (default = $w)\n\n");
		}



exit;
