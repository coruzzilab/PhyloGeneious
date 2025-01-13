#!/usr/bin/env perl

use 5.008005;
use strict;
use warnings;

use Cwd;
use Fcntl;
use AnyDBM_File;

# Local variables
my $OID_HOME;                     # OrthologID main directory
my $OID_USER_DIR;                 # OrthologID user run directory
my $OID_CONF;                     # OrthologID config file
my $OID_DATADIR;                  # OrthologID guide tree data directory
my $OID_BLASTDIR;                 # Directory to store OrthologID blast results
my $BLAST_HOME;                   # NCBI BLAST installation directory
my @INGROUP;                      # Ingroup taxa
my @OUTGROUP;                     # Outgroup taxa
my $NCPU;                         # Number of CPU to be used
my $VERBOSE;                      # Verbosity {0,1,2}
my $DEBUG = 0;                    # Debug mode
my $BLAST_RES_DB;                 # Database of all-all blast
my ( $DB_USER, $DB_PASSWORD );    # MySQL DB username and password

# BLAST e-value cutoff for potential orthologs
my $FAM_E_CUTOFF = "1e-10";

# OrthologID "global" variables
die "OID_HOME environment variable undefined!\n"     if !$ENV{'OID_HOME'};
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
my $GENE_LEN_DB  = "$OID_BLASTDIR/genelen.dbm";
my $clustersFile = "$OID_BLASTDIR/clusters";

my $unalignedFamily = "FAMILY";
my $alignedFamily   = "FAMILY.aligned";

# Initialize
#&initOID;
my $query = shift;
allBlast($query);

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
        elsif (/^\s*(MYSQL_\w+)\s*=\s*(.+?)\s*$/)
        {    # MYSQL environment (not used)
            $ENV{$1} = $2;
        }
        elsif (/^\s*DB_USER\s*=\s*(\S+)\s*$/) {    # MYSQL user (not used)
            $DB_USER = $1;
        }
        elsif (/^\s*DB_PASSWORD\s*=\s*(\S+)\s*$/) {  # MYSQL password (not used)
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
    if ( $prefix =~ /all/ ) {
        $prefix = "*";
    }

    if ( defined($prefix) && $prefix ne '*' ) {

        $blastResDB = "$OID_BLASTDIR/$prefix" . "_blastres.dbm";
        $geneLenDB  = "$OID_BLASTDIR/$prefix" . "_genelen.dbm";
    }
    else {
        print "blastall wit $prefix\n";
        $blastResDB = $BLAST_RES_DB;
        $geneLenDB  = $GENE_LEN_DB;
    }
    tie( %blastRes, "AnyDBM_File", $blastResDB, O_RDONLY, 0666 )
      or die "Cannot open $blastResDB: $!\n";

    #			tie (%geneLen, "AnyDBM_File", $geneLenDB,O_RDONLY,0666)
    #				or die "Cannot open $geneLenDB: $!\n";
    my $nk    = (%blastRes);
    my @nkeys = ();
    if ($nk) {
        @nkeys = keys %blastRes;
        $nk    = scalar(@nkeys);
        print "scalar blastres $nk\n";
        @nkeys = keys %geneLen;
        $nk    = scalar(@nkeys);
        print "scalar genelen $nk\n";
    }
    else {
        print "blastres empty\n";
        exit 5;
    }

    #			foreach my $locus ( keys %blastRes ) {
    #				$blastRes{$locus} = $blastRes{$locus};
    while ( ( my $locus, my $blastdta ) = each %blastRes ) {
        print ">$locus\n";
        print "$blastdta\n";

        #                   print "$blastRes{$locus}\n";
    }
    untie %blastRes;

    #            print "*************************\n";
    #            print "start of geneLen\n";
    #			foreach my $locus ( keys %geneLen ) {
    #				$geneLen{$locus} = $geneLen{$locus};
    #                   print ">$locus\n";
    #                   print "$geneLen{$locus}\n";
    #			}
    #			untie %geneLen;

}

