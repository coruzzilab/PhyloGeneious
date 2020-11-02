#!/usr/bin/env perl

use 5.008005;
use strict;
use warnings;

use Cwd;
use Fcntl;
use NDBM_File;

# Local variables
my $OID_HOME;        # OrthologID main directory
my $OID_USER_DIR;    # OrthologID user run directory
my $OID_CONF;        # OrthologID config file
my $OID_DATADIR;     # OrthologID guide tree data directory
my $OID_BLASTDIR;    # Directory to store OrthologID blast results
my $BLAST_HOME;      # NCBI BLAST installation directory
my @INGROUP;         # Ingroup taxa
my @OUTGROUP;        # Outgroup taxa
my $NCPU;            # Number of CPU to be used
my $VERBOSE;         # Verbosity {0,1,2}
my $DEBUG = 0;       # Debug mode
my $BLAST_RES_DB;    # Database of all-all blast
my ( $DB_USER, $DB_PASSWORD );    # MySQL DB username and password

# BLAST e-value cutoff for potential orthologs
my $FAM_E_CUTOFF     = "1e-10";


# OrthologID "global" variables
die "OID_HOME environment variable undefined!\n" if !$ENV{'OID_HOME'};
die "OID_USER_DIR environment variable undefined!\n" if !$ENV{'OID_USER_DIR'};
$OID_HOME       = $ENV{'OID_HOME'};
$OID_USER_DIR   = $ENV{'OID_USER_DIR'};
$OID_BLASTDIR   = "$OID_USER_DIR/blast";
$BLAST_RES_DB   = "$OID_BLASTDIR/blastres.dbm";
$OID_DATADIR    = "$OID_USER_DIR/data";
$ENV{'BLASTDB'} = "$OID_USER_DIR/blastdb";
$ENV{'PATH'}    = "$OID_HOME/bin:$ENV{'PATH'}";
$OID_CONF       = "$OID_USER_DIR/config";
#die "OrthologID config file not found!\n" if !-r $OID_CONF;

# Clustering
my $GENE_LEN_DB = "$OID_BLASTDIR/genelen.dbm";
my $clustersFile = "$OID_BLASTDIR/clusters";

my $unalignedFamily = "FAMILY";
my $alignedFamily = "FAMILY.aligned";

# Initialize
&initOID;
my $query = shift;
#allBlast($query);
makeTree($query);
sub initOID {

	# Defaults
	$NCPU = 1;

	# Parse configuration file
	open CONF_FH, "$OID_CONF" or die "Cannot open $OID_CONF: $!\n";
	while (<CONF_FH>) {
		chomp;
		next if (/^#/);
		if (/^\s*BLASTHOME\s*=\s*(.+?)\s*$/)
		{    # Installation directory of NCBI BLAST
			$BLAST_HOME = $1;
		}
		elsif (/^\s*INGROUP\s*=\s*(.+?)\s*$/) {    # Ingroup species prefixes
			$_       = $1;
			@INGROUP = split;
		}
		elsif (/^\s*OUTGROUP\s*=\s*(.+?)\s*$/) {    # Outgroup species prefixes
			$_        = $1;
			@OUTGROUP = split;
		}
		elsif (/^\s*NCPU\s*=\s*(\d+)/) {            # Number of CPUs available
			$NCPU = $1;
		}
		elsif (/^\s*VERBOSE\s*=\s*(\d+)/) {         # Verbosity
			$VERBOSE = $1;
		}
		elsif (/^\s*DEBUG\s*=\s*(\d+)/) {           # Debug mode
			$DEBUG = $1;
		}
		elsif (/^\s*(MYSQL_\w+)\s*=\s*(.+?)\s*$/) {    # MYSQL environment (not used)
			$ENV{$1} = $2;
		}
		elsif (/^\s*DB_USER\s*=\s*(\S+)\s*$/) {        # MYSQL user (not used)
			$DB_USER = $1;
		}
		elsif (/^\s*DB_PASSWORD\s*=\s*(\S+)\s*$/) {    # MYSQL password (not used)
			$DB_PASSWORD = $1;
		}
	}
	close CONF_FH;

	die if !( @INGROUP && @OUTGROUP );

	return 1;
}

# Preloaded methods go here.


#
# Retrieve ingroup list
#
sub getIngroup() {
	return @INGROUP;
}

#
# Retrieve outgroup list
#
sub getOutgroup() {
	return @OUTGROUP;
}

#
# Print configuration info
#
sub printConfig() {
	print "OID_HOME = $OID_HOME\n";
	print "OID_USER_DIR= $OID_USER_DIR\n";
	print "Ingroup taxa = @INGROUP\n";
	print "Outgroup taxa = @OUTGROUP\n";
}

#
# All-all blast
# Optional argument: Ingroup taxon prefix, e.g. allBlast("At")
# No argument means all ingroup taxa.
# Asterisk (*) means combine all previously generated blast results.
#
sub allBlast {
	my $prefix  = shift;
	my $verbose = $VERBOSE;
	#### E-value cutoff
	#my $eCutOff = "1e-16";
	my $eCutOff = $FAM_E_CUTOFF;
	####

	my $outputDir = $OID_BLASTDIR;    # Directory to store blast output
	my @queryGroup;                   # BLAST query taxa
	my $blastResDB;                   # DB file to store blast results
	my $geneLenDB;                    # DB file to store gene length
	my $numProc = $NCPU;
	my %blastRes;
	my %geneLen;

#    print "in blastall with $prefix\n";
    if ($prefix =~ /all/){
        $prefix = "*";
    }

	if ( defined($prefix) && $prefix ne '*' ) {
        
		$blastResDB = "$OID_BLASTDIR/$prefix" . "_blastres.dbm";
		$geneLenDB = "$OID_BLASTDIR/$prefix" . "_genelen.dbm";
    }else {
        print "blastall wit $prefix\n";
		$blastResDB = $BLAST_RES_DB;
		$geneLenDB = $GENE_LEN_DB;
	}
			tie (%blastRes, "NDBM_File", $blastResDB,O_RDONLY,0666)
			  or die "Cannot open $blastResDB: $!\n";
			tie (%geneLen, "NDBM_File", $geneLenDB,O_RDONLY,0666)
				or die "Cannot open $geneLenDB: $!\n";
            my $nk =  (%blastRes);
            if ($nk){
                print "scaler blastres $nk\n";
            }else {
                print "blastres empty\n";
                exit 5;
            }
#			foreach my $locus ( keys %blastRes ) {
#				$blastRes{$locus} = $blastRes{$locus};
            while ((my $locus,my $blastdta) = each %blastRes) {
                   print ">$locus\n";
                   print "$blastdta\n";
#                   print "$blastRes{$locus}\n";
			}
			untie %blastRes;
			foreach my $locus ( keys %geneLen ) {
#				$geneLen{$locus} = $geneLen{$locus};
                   print ">$locus\n";
                   print "$geneLen{$locus}\n";
			}
			untie %geneLen;


}


sub makeTree {

	my $dirRE = shift;
	$dirRE = '.*' if !defined($dirRE);
	my $verbose = $VERBOSE;

	# Defaults
	my $prefix    = "oid";
	#my $paup_exe  = "paup";
	my $tnt_exe = "tnt";
	my $nrat      = 10;       # number of ratchets
	my $niter     = 100;      # number of iterations per ratchet
	my $timeLimit = 10800;    # Time limit for a heuristic search
	my $ratTimeLimit = 60;    # Time limit for a search within a ratchet iteration

	my $thisFunc = ( caller(0) )[3];
	my $oldFH;
	my %seq;
	my @outgroup;           # Outgroup taxa for current family

	# Change to data dir
	my $savDir = getcwd;
	chdir $OID_DATADIR;

	# Go over each family
	my $count = 0;
	print "Making trees ...\n" if $verbose;
	foreach my $dir (<[1-9]*>) {
		next if $dir !~ /$dirRE/;
		$count++;
		chdir $dir
		  or die "$thisFunc: failed to change to family directory $dir: $!\n";
		if ( $verbose == 1 ) {
			$| = 1;    # flush STDOUT
			print "\b" x 7;
			printf "%6d", $count;
			$| = 0;
		}
		elsif ( $verbose == 2 ) {
			print "Tree search for family $dir\n";
		}

		%seq      = ();
		@outgroup = ();

		# Check nexus file or tree file does not already exist
		if ( -s $prefix . ".tre" ) {
			print "Tree file already exists ... skipping\n"
			  if $verbose;
			chdir $OID_DATADIR;
			next;
		}

		open( ALIGN, "<", $alignedFamily )
		  or die "$thisFunc: failed to open alignment: $!\n";
		my $nam;
		while (<ALIGN>) {
			chomp;
			if (/^>/) {
				$nam = substr( $_, 1 );
				$seq{$nam} = "";
			}
			elsif ( defined $seq{$nam} ) {
				$seq{$nam} .= $_;
			}
		}
		close ALIGN;
		my $length = length( $seq{$nam} ); # assume they are all the same length
		my $ntaxa  = keys %seq;
		if ($ntaxa <4){
			print "Too few taxa. TNT cannot create tree with less than 4 taxa\n"
			 if $verbose;
			chdir $OID_DATADIR;
			next;
		}
		if (-s $prefix.".nex"){ goto RUNTNT }
		open( NEX, ">", $prefix . ".nex" )
		  or die "$thisFunc: Failed to open file for writing: $!\n";
		$oldFH = select(NEX);
		print "xread\n";
		print "$length $ntaxa\n";
		print "&[proteins nogaps]\n";
		foreach ( keys %seq ) {
			if (length($_) < 32) {
				# Try to line up the characters
				print $_ . " " x ( 32 - length($_) );
			}
			else {
				print "$_ ";
			}
			print $seq{$_} . "\n";

			# Check if taxon is outgroup
			foreach my $og (@OUTGROUP) {
				if (/^$og#/) {
					push @outgroup, $_;
					last;
				}
			}
		}
		print ";\n";
		close(NEX);
		select($oldFH);
		open (PROC, ">$prefix".".proc")
		  or die "$thisFunc: Failed to open file for writing: $!\n";
		$oldFH = select(PROC);
		# TNT commands
		print "log " . $prefix . ".log;\n";
		print "taxname=;\n";
		print "taxname +64;\n";
		print "report +200\n";
		print "nstates 32;\n";
		#print "hold 1000;\n";
		if($ntaxa >10) {
			print "mxram 2000;\n";
		}
		elsif($ntaxa > 1000) {
			print "mxram 20000;\n";
		}
		print "p $prefix.nex;\n";
#		if (@outgroup) {
#			print "outgroup @outgroup;\n";
#			print "set root=outgroup;\n";
#		}
#		else {
#			print "set root=midpoint;\n" if $ntaxa > 3;
#		}
		if ( $ntaxa < 10 ) {
			print STDOUT "___full search___\n" if $verbose == 2;
			print "mult= tbr replic 20 hold 10 ;\n";
		}
		elsif ( $ntaxa < 100 ) {    # use tbr swap
			print STDOUT "___tbr sectorialSearches___\n" if $verbose == 2;
			print "mult= tbr replic 50 hold 20 ;\n";
			print "sec: xss fuse 5 drift 10 starts 10 ;\n";
		}
		else {                      # use ratchet
			print STDOUT "___ratchet___\n" if $verbose == 2;
			print "rs0;\n";
			print "col3;\n";
			print "hold 1000;\n";
			print "timeout 120:00:00;\n";
			print "rat:iter100up4do4;\n";
			print "mu=rep20ho1rat;\n";
			print "bbreak=tbr;\n";
		}
		print "taxname];\n";
		print "nel*;\n";
		print "tchoose/;\n";
		print "export - $prefix.tre;\n";
		print "log/;\n";
		print "zz/;\n";
		print "proc;\n";
		close(PROC);

		select($oldFH);

		# Execute PAUP
		
		#my $status = system("$paup_exe -n $prefix.nex 1>/dev/null 2>&1");
		#die "$thisFunc: PAUP error: $?" unless $status == 0;

		# remove temp files if any
		# Execute tnt
		RUNTNT:
#		if($ntaxa <1000){
#			my $status = system("tnt mxram 2000, p $prefix.proc");
#			die "$thisFunc: TNT error: $?" unless $status == 0;
#		}
#		else{
#			my $status = system("tnt mxram 22000, p $prefix.proc");
#			die "$thisFunc: TNT error: $?" unless $status == 0;
#		}
#		unlink <$prefix.tre.*>;
#		unlink "$prefix.tmp" if -f "$prefix.tmp";
		chdir $OID_DATADIR;
	}

	print "\n" if $verbose == 1;
	chdir $savDir;
}
