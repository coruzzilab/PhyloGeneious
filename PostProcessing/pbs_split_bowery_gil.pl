#!/usr/bin/perl

############################### This program is free software; you can redistribute it and/or modify
############################### it under the terms of the GNU General Public License as published by
############################### the Free Software Foundation; either version 2 of the License, or
############################### (at your option) any later version.
############################### 
############################### This program is distributed in the hope that it will be useful,
############################### but WITHOUT ANY WARRANTY; without even the implied warranty of
############################### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
############################### GNU General Public License for more details.
############################### 
############################### You should have received a copy of the GNU General Public License
############################### along with this program; if not, write to the Free Software
############################### Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
###############################
###############################
###############################
############################### Copyright 2013 Damon P. Little
############################### Adapted to cluter submission by splitting PBS calculations across N
############################### jobs where N is the number of internal nodes in the tree.
############################### 2014 Kranthi Varala. Copyright still belongs to Damon P. Little
############################### Adjusted checks and fixed node matching.
############################### 2024 Veronica Sondervan.
my $l = (); #logfile
my $m = (); #matrix file
my %n = (); #nel tree file
my $t = (); #tre tree file
my $skip = 0;
for(my $k = $#ARGV; $k >= 0; $k--){
	if($ARGV[$k] eq '-l'){
		if(-e $ARGV[$k+1]){
			print(STDERR "Logfile $l will be overwritten!\n");
		}
		if(length($ARGV[$k+1])){
			$l = $ARGV[$k+1];
		}
		next;	
	}
	if($ARGV[$k] eq '-m'){
		if(-e $ARGV[$k+1]){
			$m = $ARGV[$k+1];
		}
		next;	
	}
	if($ARGV[$k] eq '-n'){
		if(-e $ARGV[$k+1]){
			$n = $ARGV[$k+1];
		}
		next;	
	}
	if($ARGV[$k] eq '-t'){
		if(-e $ARGV[$k+1]){
			$t = $ARGV[$k+1];
		}
		next;	
	}
	if($ARGV[$k] eq '-s'){
		$skip=1;
		next;	
	}
}

############################### MODULE
use Text::Balanced qw(extract_bracketed);
use Cwd;

my $OID_USER_DIR;
our $HPC='P';
die "OID_USER_DIR environment variable undefined!\n" if !$ENV{'OID_USER_DIR'};
$OID_USER_DIR   = $ENV{'OID_USER_DIR'};

# Parse configuration file
open CONF_FH, "$OID_USER_DIR/config" or die "Cannot open $OID_USER_DIR/config: $!\n";
while (<CONF_FH>) {
    chomp;
    next if (/^#/);
    if (/^\s*HPC\s*=\s*([PSC]).*$/) 
    {
        $HPC = $1;
    }
}

if(length($l) && length($m) && length($n) && length($t)){

	############################### PROC FILE
	my $p = 0;
	open(INFILE, $m) || die("Could not open $m!");
	while(my $line = <INFILE>){
		chomp($line);
		if(($line =~ m/^&\[/) && ($line =~ m/\]$/)){
			$p++;
		}
	}
	close(INFILE);
	
	my %clades = ();
	my $terminals = 0;
	open(INFILE, $n) || die("Could not open $n!");
	while(my $line = <INFILE>){
		chomp($line);
		if(length($line) && ($line !~ m/^tread/) && ($line ne 'proc-;')){
			$line =~ tr/0123456789 \(\)//cd;
			my $remainder = substr($line, index($line, '('));
			$terminals = ($line =~ tr/ / /);
			while(($remainder =~ m/\(/) && ($remainder =~ m/\)/)){
				my @bits = extract_bracketed($remainder, '()');
				$bits[0] = substr($bits[0], 1, (length($bits[0])-2));
				if(($bits[0] =~ m/\(/) && ($bits[0] =~ m/\)/)){
					$remainder = substr($bits[0], index($bits[0], '(')) . $bits[1];
					} else {
						$remainder = $bits[1];
						}
				$bits[0] =~ tr/0123456789 /0123456789 /cds;
				$bits[0] =~ s/ $//g;
				$clades{$bits[0]} = scalar(split(/ /, $bits[0]));				
			}
			last;
		}
	}
	close(INFILE);
	my @clades = keys %clades; #keys order is inconsistent
#	print join(", ", @clades)."\n";
#	print "$clades{'16 19'}\n"; #index is consistent
	#my @jobs;
	my $pwd = getcwd();
	my $tnt = "log unconstrained.log;\n";
	$tnt .= "mxram 32000;\n"; #changed to more reasonable amount
	$tnt .= "p $m;\n";
	$tnt .= "p $t;\n";
	$tnt .= "silent=buffer;\n";
	for(my $k = $p; $k > 0; $k--){
		$tnt .= "quote UNCONSTRAINED PARTITON $k;\n";
		$tnt .= "cc].;\n";
		$tnt .= "cc[\@$k.;\n";
		$tnt .= "len;\n";
	}
	$tnt .= "silent-buffer;\n";
	$tnt .= "log/;\n";
	$tnt .= "zz;\n";
	$tnt .= "proc/;\n";
	unless(-e "unconstrained.proc"){
		open(OUTFILE, ">unconstrained.proc") || die("Could not open unconstrained.proc!");
		print(OUTFILE $tnt);
		close(OUTFILE);
		open(JOBFILE,">unconstrained.sbatch") || die("Could not open unconstrained.sbatch!");
		print JOBFILE "#!/bin/bash\n#SBATCH --verbose\n#SBATCH --export=ALL\n#SBATCH --job-name=PBS.$c\n#SBATCH --output=PBS.$c.out\n#SBATCH --error=PBS.$c.err\n#SBATCH --time=24:00:00\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=32GB\n";#reduced time from 168hrs
		print JOBFILE "cd $pwd\n";
		print JOBFILE "if ! [[ -f \$HOME/.passwordfile.tnt ]]; then cp \$OID_HOME/.passwordfile.tnt \$HOME; fi\n";
		print JOBFILE "tnt bground p unconstrained.proc\n";
		close JOBFILE;
		if($skip < 1){
			if($HPC=="S"){
				system "sbatch unconstrained.sbatch";
			}
			else{
				system "qsub unconstrained.sbatch";
			}
		}
	}
	for(my $c=$#clades;$c >=0; $c--){
#		if($clades{$clades[$c]} < ($terminals-1)){
			my $tnt = "log $l.$c;\n";
			$tnt .= "mxram 32000;\n"; #changed to more reasonable amount
			$tnt .= "p $m;\n";
			$tnt .= "p $t;\n";
			$tnt .= "silent=buffer;\n";
			$tnt .= "quote BEGIN PBS;\n";
			$tnt .= "hold 1000;\n";
			$tnt .= "rat:iter50up4do4;\n";
			$tnt .= "mult:rep1ho1rat;\n";
			$tnt .= "k0;\n";
			$tnt .= "quote CLADE TO CONSTRAIN IS $clades[$c];\n";
			$tnt .= "force-[$clades[$c]];\n";
			$tnt .= "const=;\n";
			$tnt .= "cc[.;\n";
			$tnt .= "mult;\n";
			$tnt .= "bb;\n";
			for(my $j = $p; $j > 0; $j--){
				$tnt .= "quote CONSTRAINED PARTITON $j NODE $c;\n";
				$tnt .= "cc].;\n";
				$tnt .= "cc[\@$j.;\n";
				$tnt .= "len;\n";
			}
			$tnt .= "quote END PBS;\n";
			$tnt .= "silent-buffer;\n";
			$tnt .= "log/;\n";
			$tnt .= "zz;\n";
			$tnt .= "proc/;\n";
			unless(-e "temporary-tnt.$c.proc"){
				open(OUTFILE, ">temporary-tnt.$c.proc") || die("Could not open temporary-tnt.$c.proc!");
				print(OUTFILE $tnt);
				close(OUTFILE);
				open JOBFILE, ">pbs.$c.sub" || die "Could not open pbs.$c.sub";
				print JOBFILE "#!/bin/bash\n#SBATCH --verbose\n#SBATCH --export=ALL\n#SBATCH --job-name=PBS.$c\n#SBATCH --output=PBS.$c.out\n#SBATCH --error=PBS.$c.err\n#SBATCH --time=24:00:00\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=32GB\n";#reduced time from 168hrs
				print JOBFILE "cd $pwd\n";
				print JOBFILE "if ! [[ -f $HOME/.passwordfile.tnt ]]; then cp $OID_HOME/.passwordfile.tnt $HOME; fi\n";
				print JOBFILE "tnt bground p temporary-tnt.$c.proc\n";
				close JOBFILE;
			}
			if($skip < 1){
				if($HPC=="S"){
					my $jid = `sbatch pbs.$c.sub`;
				}
				else{
					my $jid = `qsub pbs.$c.sub`;
				}
				#push @jobs, $jid;
				sleep(90); #prevent TNT from overloading
			}
#		}
	}

	############################### TNT
	#qx/tnt bground p temporary-tnt.proc/;
	#unlink(temporary-tnt.proc);
	###############################
	
	############################### Check Job completion
	while($skip < 1){
		if($HPC=="S"){
			my $jcnt = `squeue -u $USER |grep PBS|wc -l`;
		}
		else{
			my $jcnt = `qstat -u $MYUSER |grep PBS|wc -l`;
		}
		chomp $jcnt;
		if($jcnt<1){last;}
		sleep(60);
	}

	###############################
	system "cat unconstrained.log >$l";
	system "cat $l.[0-9]* >>$l";
	############################### READ LOGS
	my $buffer = "clade,partition,pbs\n";
	my @numbertonode = ();
	$#numbertonode=$#clades;
	for(my $c=$#clades;$c >=0; $c--){
		open(INFILE, "temporary-tnt.$c.proc") || die("Could not open temporary-tnt.$c.proc !");
		while(my $line = <INFILE>){
			chomp($line);
			$line =~ tr/ / /s;
			$line =~ s/^ //g;
			if(length($line)){
				if($line =~ /force/){#-\[(.+)\];/){
					$line =~ s/force-\[(.+)\];/\1/;
					for(my $ck=$#clades;$ck >=0; $ck--){
						if("$line" eq "$clades[$ck]"){$numbertonode[$c] = $ck;}
					}
#					print "$c | $line | $numbertonode[$c] | $clades[$numbertonode[$c]]\n";
				}
			}
		}
		close(INFILE);
	}
#	print join(", ", @numbertonode)."\n";
	open(INFILE, $l) || die("Could not open $l !");
	my $constrained = -1;
	my $k = 0;
	my @lengths = ();
	my $node = ();
	my @partitions = ();
	my $tree = 0;
	my $unconstrained = -1;
	while(my $line = <INFILE>){
		chomp($line);
		$line =~ tr/ / /s;
		$line =~ s/^ //g;
		if(length($line)){
			my @bits = split(/ /, $line);
			if($tree && ($k == 2)){
				push(@lengths, @bits[1..$#bits]);
				next;
			}
			if($tree && ($k == 3)){
				$tree = 0;
				my $mean = ();
				for(my $j = $#lengths; $j >= 0; $j--){
					$mean += $lengths[$j];
				}
				$mean = $mean / ($#lengths+1);
				if(($unconstrained != -1) && ($mean >= 0)){
					$partitions[$unconstrained] = $mean;
				}
				elsif(($constrained != -1) && ($mean >= 0)){
						$buffer .= $clades[$node] . ',' . $constrained . ',' . sprintf('%.4f', ($mean - $partitions[$constrained])) . "\n";
						if(length($buffer) > 10000){
							print($buffer);
							$buffer = ();
						}
				}
				else {
						print(STDERR "Error reading TNT log file!\n");
						last;
				}
			}
			if(($bits[0] eq 'UNCONSTRAINED') && ($bits[1] eq 'PARTITON')){
				$unconstrained = $bits[2];
				next;
			}
			if(($bits[0] eq 'CONSTRAINED') && ($bits[1] eq 'PARTITON')){
				$constrained = $bits[2];
				$node = $numbertonode[$bits[4]];
				$unconstrained = -1;
				next;
			}
			if(($bits[0] eq 'Tree') && ($bits[1] eq 'lengths')){
				$k = 0;
				$tree = 1;
				@lengths = ();
				next;
			}
		}
		else {
			$k++;
		}
	}
	close(INFILE);
	print($buffer);
}
else { ############################### PRINT USAGE INFO
	print("\nA PERL script for computing Partition Bremer Support with TNT.\n");
	print("USAGE: pbs.pl -l logfile.log -m matrix.tnt -n consensus.nel -t trees.tre [-s skip job submission]\n\n");
	print("WHERE: Matrix is a TNT formatted matrix with one partion per data block.\n");
	print("The consensus tree and the most parsimonious trees must be in numeric format.\n\n");
}

