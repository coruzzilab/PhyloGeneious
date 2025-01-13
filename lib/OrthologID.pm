#
# OrthologID Perl module
#
# Copyright (C) 2006-2012 Ernest K. Lee
#
#    This file is part of OrthologID.
#
#    OrthologID is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    OrthologID is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with OrthologID.  If not, see <http://www.gnu.org/licenses/>.
#
# NOTE:
# - The OID_HOME and OID_USER_DIR environment variables must be defined.
# - Gene names must be prefixed by species names as defined in config, in the format
#   species#geneID, e.g. Arath#At1g12345 or SpX#NP_123456.1 (DO NOT USE "-" in species/gene names).
# - The command line 'paup' and 'mcl' executables must be in your PATH.
#
# Author: Ernest K Lee <elee@amnh.org>
#

package OrthologID;

use 5.008005;
use strict;
use warnings;

use Cwd;
use Fcntl;
use Time::Piece;
use AnyDBM_File;
use DB_File;

#use IO::Compress::Gzip qw(gzip $GzipError);
use OrthologID::PhyloTree;
use OrthologID::Utils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
    'all' => [
        qw(

        )
    ]
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
  printConfig
  getIngroup
  getOutgroup
  allBlast
  makeFamily
  rdgeners
  nextblst
  mclCluster
  mclmcx
  make1fam
  chkuptm
  setglsttm
  alignFamily
  makeTree
  schedatqs
  findOrthologs
  %procptr
  %limval
  $TNTA
  fmthpc
  getqids
  fixqs
  testqs
);

#use vars qw(%procptr $TNTA);
our $VERSION = '0.01';

# Local variables
my $OID_HOME;                     # OrthologID main directory
my $OID_USER_DIR;                 # OrthologID user run directory
my $OID_WRAPPER = "";                  # OrthologID container wrapper script
my $OID_CONF;                     # OrthologID config file
my $OID_DATADIR;                  # OrthologID guide tree data directory
our $OID_BLASTDIR;                # Directory to store OrthologID blast results
my $BLAST_HOME;                   # NCBI BLAST installation directory
our $SEARCHTYPE = 'B';            # similarity search method default B=BLASTP,(D=DIAMOND, M=mmseqs2)
our $TREEPROGRAM = 'TNT';         # tree building program default TNT,(OBLONG)
our $OID_COMB;                    # combined fa files for all input
my $OID_PROC;                     # oid.proc input file
our @INGROUP;                     # Ingroup taxa
our @OUTGROUP;                    # Outgroup taxa
our $NCPU;                        # Number of CPU to be used
our $VERBOSE;                     # Verbosity {0,1,2}
my $DEBUG = 0;                    # Debug mode
our $HPC  = 'P';    #kind of hpc system default P=PBS,(S=SLUrm, C=CGE)
our $TNTA = 200;    # cutoff for serial tnt processing
our $TNTB = 500;    # cutoff for intermediate (multi-thread) tnt processing
our $myjob;
my $TNTB_THR = 10;      # number of threads to run intermediate tnt
my $TNTC_THR = 10;      # number of singleton threads (used to be maxqs)
my $TNTC     = 1000;    #big singleton
our $TNTD = 10000;      #print but don't bother - I think it cannot be done
my $BLAST_RES_DB;       # Database of all-all blast
our %limval   = ();     #limval{cutoff} = [gb, hours]
our @limkey   = ();     #sorted keys of limval
our %sztaxa   = ();
our $BLST_MIN = 60;     ## this is the avg time of a blast task in min
my ( $DB_USER, $DB_PASSWORD );    # MySQL DB username and password
our $MAXQS = $TNTC_THR;           ## make global
my $MAXFORK = 2;                  ##make global
my $MYUSER  = `whoami`;
our $glsttm;
my $printct  = 0;
my $SCHEDONE = 'log/job/schedone';
my $SCHEDLOG = 'log/job/schedlog.txt';
my $BMDONE   = 'log/job/bigmdone';
my $PLDONE   = 'log/job/pooldone';
#our $JOBCAP = 2000; redundant with MAXQS

# BLAST e-value cutoff for potential orthologs
my $FAM_E_CUTOFF = "1e-10";

# OrthologID "global" variables
die "OID_HOME environment variable undefined!\n"     if !$ENV{'OID_HOME'};
die "OID_USER_DIR environment variable undefined!\n" if !$ENV{'OID_USER_DIR'};
$OID_HOME       = $ENV{'OID_HOME'};
$OID_USER_DIR   = $ENV{'OID_USER_DIR'};
if ($ENV{'OID_WRAPPER'}) {
    $OID_WRAPPER    = $ENV{'OID_WRAPPER'};
}
$OID_BLASTDIR   = "$OID_USER_DIR/blast";
$BLAST_RES_DB   = "$OID_BLASTDIR/blastres.blst";
$OID_DATADIR    = "$OID_USER_DIR/data";
$ENV{'BLASTDB'} = "$OID_USER_DIR/blastdb";
$ENV{'PATH'}    = "$OID_HOME/bin:$ENV{'PATH'}";
$OID_CONF       = "$OID_USER_DIR/config";
$OID_PROC       = "$OID_USER_DIR/procfiles.txt";
$OID_COMB       = "$OID_USER_DIR/blastdb/combined.fa";

die "OrthologID config file not found!\n" if !-r $OID_CONF;

# Clustering
my $GENE_LEN_DB  = "$OID_BLASTDIR/genelen.blst";
my $clustersFile = "$OID_BLASTDIR/clusters";
my $clusterdone  = "$OID_BLASTDIR/.mci.done";
my $mcxdone      = "$OID_BLASTDIR/.mcx.done";
my $familydone   = "$OID_DATADIR/.family.done";

my $unalignedFamily = "FAMILY";
our $alignedFamily = "FAMILY.aligned";
## Proc files for tnt
our %procptr = ()
  ; #entry is cutoff, points to array with proc or hash ncpu with array of arrays
our @cutoff = ();    #sorted cutoffs
our @cutmax = ();    #number of procs ie if multitask need 1 proc for each task

# Initialize
&initOID;

#push @EXPORT,qw (%procptr);
#push @EXPORT,qw ($TNTA);
sub init_myjob {

    if ( $HPC =~ /S/ ) {
        $myjob = $ENV{"SLURM_JOB_ID"};
        print "my job $myjob\n" if ($myjob);
    }
    else {
        $myjob = $ENV{"PBS_JOBID"};
        print "my job $myjob\n" if ($myjob);
    }
}

sub init_proc {
    my $readproc = 0;
    my $curproc  = 0;
    my $maxproc  = 0;
    my $curcut   = 0;
    my $ar       = 0;
    open PROC, "$OID_PROC" or die "Cannot open $OID_PROC: $!\n";
    while (<PROC>) {
        chomp;
        if ($readproc) {
            push @$ar, $_;
            if (/proc\;/) {
                $readproc = 0;

          #                print "cutoff $curcut, proc $curproc max $maxproc\n";
          #                print join "\n",@$ar;
            }
        }
        else {
            if (/\s*CUTOFF\s*=\s*(\d+)/) {
                $curcut = $1;
                if ( exists $procptr{$curcut} ) {
                    print "duplicate cutoff $_\n";
                    next;
                }
                $curproc  = 0;
                $maxproc  = 1;
                $readproc = 1;
                if (/\s*NPROC\s*=(\d*)/) {
                    $maxproc = $1;
                }
                $procptr{$curcut} = [$maxproc];
                $curproc++;
                push @{ $procptr{$curcut} }, [];

                #                $ar = \{${$procptr{$curcut}}[$curproc];
                $ar = $procptr{$curcut}->[$curproc];
            }

            elsif (/\s*PROC\s*/) {
                $curproc++;
                if ( $curproc > $maxproc ) {
                    print "curproc>maxproc $maxproc for cutoff $curcut\n";
                    next;
                }
                push @{ $procptr{$curcut} }, [];
                $ar       = $procptr{$curcut}->[$curproc];
                $readproc = 1;
            }
        }
    }    #read file
    @cutoff = sort { $a <=> $b } keys %procptr;
    for ( my $i = 0 ; $i <= $#cutoff ; $i++ ) {
        $cutmax[$i] = ${ $procptr{ $cutoff[$i] } }[0];
        print "cutoff $cutoff[$i], max $cutmax[$i]\n";
    }
}

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
        elsif (/^\s*HPC\s*=\s*([PSC]).*$/) {
            $HPC = $1;
        }
        elsif (/^\s*SEARCHTYPE\s*=\s*([BDM]).*$/) {   # search method
            $SEARCHTYPE = $1;
        }
        elsif (/^\s*TREEPROGRAM\s*=\s*(.+?)\s*$/) {   # tree build method
            $TREEPROGRAM = $1;
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
        elsif (/^\s*BLSTMIN\s*=\s*(\d+)/) { #avg time in minute to set for blast
            $BLST_MIN = $1;
        }
        elsif (/^\s*TNTA\s*=\s*(\d+)/) {    #serial tnt processing size cutoff
            $TNTA = $1;
        }
        elsif (/^\s*TNTB\s*=\s*(.+?)\s*$/)
        {    #intermediate tnt group cutoff and nthreads
            $_ = $1;
            my @tmp = m/(\d+)/g;
            $TNTB     = $tmp[0];
            $TNTB_THR = $tmp[1];
        }
        elsif (/^\s*TNTC_THR\s*=\s*(\d+)/) {    #singlton tnt max jobs
            $TNTC_THR = $1;
            $MAXQS    = $TNTC_THR;
        }
        elsif (/^\s*TNTD\s*=\s*(\d+)/) {    #too big tnt processing size cutoff
            $TNTD = $1;
        }
        elsif (/^\s*LIMIT\s*=\s*(.+?)\s*$/) {    #singlton tnt max jobs
            $_ = $1;
            my @tmp = split;
            $limval{ $tmp[0] } = [ $tmp[1] ];
            push @{ $limval{ $tmp[0] } }, $tmp[2];
            print "entered LIMIT ", join ",", @tmp, "\n";

        }

        elsif (/^\s*(MYSQL_\w+)\s*=\s*(.+?)\s*$/)
        {                                        # MYSQL environment (not used)
            $ENV{$1} = $2;
        }
        elsif (/^\s*DB_USER\s*=\s*(\S+)\s*$/) {    # MYSQL user (not used)
            $DB_USER = $1;
        }
        elsif (/^\s*DB_PASSWORD\s*=\s*(\S+)\s*$/) {  # MYSQL password (not used)
            $DB_PASSWORD = $1;
        }
        elsif (/^\s*MAXQS\s*=\s*(\d+)/) {            # Maximum qsub jobs
            $MAXQS    = $1;
            $TNTC_THR = $MAXQS;                      #synonym
        }
    }
    close CONF_FH;
    @limkey = keys %limval;

    #    print @limkey," limkey\n";
    print "maxqs $MAXQS ncpu $NCPU\n";
    if ( !$#limkey ) {    #empty
        $limkey[0] = "$TNTA";
        $limval{$TNTA} = [ 12, 12 ];    #12gb, 12 hours
        print "no limit entered useing 200,12,12 (12 gb,12 hours)\n";
    }
    else {
        @limkey = sort { $a <=> $b } @limkey;
        foreach my $v (@limkey) {
            my $str1 = $limval{$v}[0];
            my $str2 = $limval{$v}[1];
            print "limit $v is  $str1", "gb , $str2 hours\n";
        }
    }

    die if !( @INGROUP && @OUTGROUP );
    init_proc();     #read in proc files
    init_myjob();    #setup myjob

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

sub fixqs {
    my $qs = shift;
    if ( $qs !~ /^\d+$/ ) {
        my @stf = split /\s/, $qs;
        if ( $stf[-1] =~ /^\d+$/ ) {
            $qs = $stf[-1];
            return $qs;
        }
        else {
            print "truble breaking up $qs\n";
            return 0;
        }
    }
    return $qs;
}    #end fixqs

sub fmthpc {
    ( my $gb, my $wall, my $cpus, my $tasks, my $script, my @args ) = @_;
    if ( $HPC =~ /S/ ) {
        my $PARAMS =
          "-N 1 --ntasks-per-node $tasks -c $cpus -t $wall --mem $gb ";
        my $argstr = '';
        foreach my $arg (@args) {
            $argstr .= " $arg";
        }
        my $submit = q/sbatch / . "$PARAMS $OID_WRAPPER $script $argstr";
        return $submit;
    }
    else {    #HPC = 'P' or undefined
        my $walmem = "mem=$gb" . ",walltime=" . "$wall";
        my $PARAMS = "-l nodes=1:ppn=$NCPU -l $walmem";
        my $argstr = "";
        my $narg   = @args;
        if ( @args > 0 ) {
            $argstr = "-v ";
            my $argnum = 0;
            foreach my $arg (@args) {
                $argnum++;
                $argstr .= "arg$argnum" . "=$arg";
                $argstr .= "," if ( $argnum < $narg );
            }
        }
        my $submit = q/qsub  / . "$PARAMS $OID_WRAPPER $script $argstr";
        return ($submit);
    }
}

sub getqids {
    my $rqs   = shift;
    my $flg   = 1;
    my @stuff = ();
    my $qs;

    if ( $HPC =~ /S/ ) {
        my $pid = open( QS, "squeue -l -u $MYUSER |" );
        my $jnk = <QS>;
        $jnk = <QS>;
        while (<QS>) {
            chomp;
            @stuff = split;
            if ($flg) {

                #                print " squeue " ,join " ",@stuff;
                #                print "squeue ",scalar(@stuff);
                #                print "\nX$stuff[0]X\n";
                #                $flg = 0;
            }
            next if ( $stuff[0] =~ /\D/ );
            $qs = "$stuff[0]";
            $$rqs{"$qs"} = 1;             #mark active
            if ( $printct < 5 ) {
                print "squeue ", scalar(@stuff);
                print "  $qs and $$rqs{$qs}\n";
            }
            $printct++;
            $printct = 0 if $printct > 100;
        }

        return;
    }
    else {

        my $pid = open( QS, "qstat -e -u $MYUSER |" );
        my $jnk = <QS>;
        $jnk = <QS>;

        #    print "testqs $jnk\n";
        while (<QS>) {
            chomp;
            @stuff = split;
            next if ( $stuff[0] =~ /\D/ );
            if ( $stuff[-2] =~ /R|Q/ ) {
                $qs = "$stuff[0]";
                $$rqs{"$qs"} = 1;             #mark active
            }
        }
    }
}

sub setglsttm {
    my $t = time();
    $glsttm = $t;
    if ( -s 'starttime' ) {
        open FI, "<starttime";
        my $x = <FI>;
        close FI;
        if ($x) {
            chomp($x);
            my $dif = $t - $x;
            my $xx  = localtime($x);
            print "$xx in starttime used to time us, dif $dif \n";
            $glsttm = $x;
        }
    }
}

sub chkuptm {
    my $tmnw = time();
    my $td   = $tmnw - $glsttm;
    if ( $td > 72000 ) {
        print "time now ", localtime($tmnw), " diff ", $td, "\n";
        restartme();    #restarts my task
    }
}

#
# All-all blast
# Optional argument: Ingroup taxon prefix, e.g. allBlast("At")
# No argument means all ingroup taxa.
# Asterisk (*) means combine all previously generated blast results.
#
# step 1 just convert from blastall to blastp (blast+)
#  it is much more efficient and logical to run entire (combimed) input again each database
#  rather than run each fasta file against all databases - especially if database very large
#  because then db cashed and read in only 1x
#  -- this is just proof of concept to see if blastp works
#  step 2 - change order of blast: combined.fa has all fasta records
#  makes sense to run this file again 1 db at a time (so blastp loaded once, db is cached)
#  so $prefix will have db i run everything against
sub allBlast {
    my $prefix = shift;
    our $verbose = $VERBOSE;
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

    if ( !-d $outputDir ) {
        mkdir $outputDir, 0755 or die "Cannot create directory: $!\n";
    }
    if ( defined($prefix) && $prefix ne '*' ) {

        #        $blastResDB = "$OID_BLASTDIR/Part$prefix" . "_blastres.dbm";
        #        $geneLenDB = "$OID_BLASTDIR/Part$prefix" . "_genelen.dbm";
        @queryGroup = ($prefix);
    }
    else {
# here we normally combine blast databases into big blastresdb.blst, also create genLendb
# problem - file may be very very big - and die on memory  or just be slow
# actual mclCluster just uses 1st target of each gene - so we don't need to write
# all the rest - makes for much smaller file and array in memory
        $blastResDB = $BLAST_RES_DB;
        $geneLenDB  = $GENE_LEN_DB;

        #@queryGroup = (@INGROUP, @OUTGROUP);
        if ( $prefix eq '*' ) {
            my $Parts = `ls $OID_BLASTDIR/Part*faa | wc -l`;
            chomp $Parts;

            # Combine previously generated blast results for all taxa
            warn "$blastResDB exists.  Previous results will be overwritten.\n"
              if -f $blastResDB;

           #            tie (%blastRes, "AnyDBM_File", $blastResDB,O_WRONLY|O_CREAT,0666)
           #              or die "Cannot open $blastResDB: $!\n";
           #            tie (%geneLen, "AnyDBM_File", $geneLenDB,O_WRONLY|O_CREAT,0666)
           #                or die "Cannot open $geneLenDB: $!\n";
           # this part creates geneLenDB - ok don't change
            open FH, "$OID_COMB"   or die "cannot open combdb.fa\n";
            open FA, ">$geneLenDB" or die "cannot open genelendb $geneLenDB\n";
            my $geneid = "";
            my $genSeq = "";
            while (<FH>) {
                if (/^>(\S+)/) {

       #                   $geneLen{$geneid} = length($genSeq) if ($geneid);
       #                   print FA,"$geneid,\t,length($genSeq)\n" if ($geneid);
                    if ($geneid) {
                        my $genlen = length($genSeq);
                        print FA "$geneid\t$genlen\n";
                    }
                    $geneid = $1;
                    $genSeq = "";
                }
                else {
                    $genSeq .= $_;
                }
            }    #while fh
            if ($geneid) {
                my $genlen = length($genSeq);
                print FA "$geneid\t$genlen\n";
            }
            close FA;

            # now for change as described above - much more interleaved
            open FA, ">$blastResDB"
              or die "cannot open blastresdb $blastResDB\n";

            #            $geneLen{$geneid} = length($genSeq) if ($geneid);
            #            my %allRes;
            #            my %allLen;

            my $locus  = '';
            my $target = '';
            my $count  = 0;
            my $ct     = 0;

            #        foreach my $sp (@INGROUP, @OUTGROUP) {
            for ( my $part = 1 ; $part <= $Parts ; $part++ )
            {    # each part has a unique set of source genes
                $ct       = localtime;
                %blastRes = ();
                print "$ct Combining BLAST and length results for $part\n"
                  if $verbose;
                $| = 1;    # flush STDOUT

        #               for(my $part=1;$part<=$Parts;$part++){ #$sp (@INGROUP, @OUTGROUP)
                foreach my $sp ( @INGROUP, @OUTGROUP ) {

            #                my %spGeneLen;
            #                my $spBlastResDB = "$OID_BLASTDIR/Part$part" . "_blastres.dbm";
                    my $spblstfl = "$OID_BLASTDIR/$sp.Part$part.blst";
                    if ( !-f $spblstfl ) {
                        print STDERR
                          "Blast results for $sp.Part$part not found!\n";
                        next;
                    }
                    open FH, "$spblstfl"
                      or print "blast failed for $spblstfl\n";
                    my $oldlocus = '';
                    my %tseen    = ();
                    while (<FH>) {
                        my @stuff = split;
                        $target = $stuff[1];    #target gene
                        $locus  = $stuff[0];
                        if ( $locus ne $oldlocus ) {
                            %tseen            = ();
                            $oldlocus         = $locus;
                            $blastRes{$locus} = ""
                              if !defined $blastRes{$locus};
                        }
                        next if defined $tseen{$target};
                        $blastRes{$locus} .= $_;
                        $tseen{$target} = 1;

                        #                    $spBlastRes{$locus} .= $_;
                    }
                    close FH;
                }
                foreach $locus ( keys %blastRes ) {
                    print FA ">$locus\n$blastRes{$locus}"
                      ;    #i think there are imbeded c/r
                    $count++;
                }
            }
            close FA;
            $ct = localtime;
            print "$ct completed blastres merge for $count files\n";

            #                my $spGeneLenDB = "$OID_BLASTDIR/Part$part" . "_genelen.dbm";
            #                if ( !-f $spGeneLenDB ) {
            #                    print STDERR "Gene lengths for Part$part not found!\n";
            #                    next;
            #                }

#                print "Combining BLAST and length results for $part\n";
#                tie (%spBlastRes, "AnyDBM_File", $spBlastResDB,O_RDONLY,0666);
#                foreach my $locus ( keys %spBlastRes ) {
#                    $allRes{$locus} = $spBlastRes{$locus};
#                }
#                untie %spBlastRes;
#                tie (%spGeneLen, "AnyDBM_File", $spGeneLenDB,O_RDONLY,0666);
#                foreach my $locus ( keys %spGeneLen ) {
#                    $allLen{$locus} = $spGeneLen{$locus};
#                }
#                untie %spGeneLen;
#            print "All BLAST and length results in memory. Now writing to disk\n" if $verbose;
#            foreach my $locus(keys %allRes){
##                $blastRes{$locus} = $allRes{$locus};
            #                $geneLen{$locus} = $allLen{$locus};
            #            }
            #            untie %blastRes;
            #            untie %geneLen;
            return 1;
        }
    }

    my ( $BLASTALL, $FASTACMD );
    $BLASTALL = "blastp";

#    $BLASTALL = defined($BLAST_HOME) ? "$BLAST_HOME/bin/blastall" : "blastall";
#    $FASTACMD = defined($BLAST_HOME) ? "$BLAST_HOME/bin/fastacmd" : "fastacmd";
# Blast each query sequence against everything
# the following copying to temp dir I think is unnessary and actually takes time
#    my $tmpDir;
#    if (defined $ENV{'PBS_JOBTMP'}) {$tmpDir = $ENV{'PBS_JOBTMP'}."/oid-$$";}
#    else{$tmpDir= "/tmp/oid-$$";}
#    if ( !-d $tmpDir ) {
#        mkdir $tmpDir, 0755 or die "Cannot create $tmpDir: $!";
#    }
#    my $tmpFile = "$tmpDir/Part$prefix.faa";
#    system "cp $OID_BLASTDIR/Part$prefix.faa $tmpFile";
## end of copy - commented out
    my $blstinfile = "$OID_BLASTDIR/Part$prefix.faa";    #instead use original
    my ( $blastCmd, @blastOut );

    #    my $count;

    #    warn "$blastResDB exists.  Previous results will be overwritten.\n"
    #      if -f $blastResDB;
    #    tie (%blastRes, "AnyDBM_File", $blastResDB,O_WRONLY|O_CREAT,0666)
    #        or die "Cannot open $blastResDB: $!\n";
    #    tie %blastRes, "AnyDBM_File", $blastResDB
    #      or die "Cannot open $blastResDB: $!\n";
    #    tie (%geneLen, "AnyDBM_File", $geneLenDB,O_WRONLY|O_CREAT,0666)
    #      or die "Cannot open $geneLenDB: $!\n";
    #    open FA, "$tmpFile";
    #    open DBG, ">$OID_BLASTDIR/checkLength.$prefix.txt";
##    my ($SeqID,$Seq);
    #    $SeqID = '';
    #    while(my $l = <FA>){
    #        chomp $l;
    #         if ($l =~ /^>(\S+)/) {
    #            if ($SeqID) {
    #                  $geneLen{$SeqID}=length($Seq);
    #                print DBG "$SeqID\t$geneLen{$SeqID}\n";
    #            }
    #            $SeqID=$1;
    #            $Seq="";
    #        }
##        else{$Seq.=$l;}
    #    }
    #    $geneLen{$SeqID}=length($Seq);
    #    print DBG "$SeqID\t$geneLen{$SeqID}\n";
    #    close DBG;
    #    close FA;
    print "BLASTing .. $prefix\n";
    foreach my $sp ( @INGROUP, @OUTGROUP ) {
        my $t = localtime;
        print "$t db sb $sp\n";
        my $outfile = "$OID_BLASTDIR/$sp.Part$prefix.blst";
        if ( system( "grep -q '\\.$sp\\.Part$prefix\\.done' $OID_BLASTDIR/.$sp.Part.done") < 1 ) {
            print "$sp.Part$prefix done - skipping\n" if $verbose > 1;
            next;
        }
# Run blast
#        $blastCmd =
#            "$BLASTALL -a $numProc -p blastp -d $sp -I -e $eCutOff -m 8 -i $tmpFile";
            if ($SEARCHTYPE eq 'B') {
                $blastCmd =
                    "$BLASTALL -db $sp -query $blstinfile -outfmt 6 -evalue $eCutOff "
                    . "-num_threads $numProc  -out $outfile";
            }
            elsif ($SEARCHTYPE eq 'D') {
                my $spdb = $OID_USER_DIR . "/blastdb/" . $sp . ".dmnd";
                # $eCutOff = 
                $blastCmd =
                    "/scratch/vms351/diamond $BLASTALL --db $spdb --query $blstinfile --evalue $eCutOff "
                    . "--threads $numProc --out $outfile --outfmt 6 --ultra-sensitive --max-target-seqs 500";
            }
            elsif ($SEARCHTYPE eq 'M') {
                my $spdb = $OID_USER_DIR . "/blastdb/" . $sp . ".DB";
                # $eCutOff = 
                $blastCmd =
                    "/scratch/vms351/mmseqs/bin/mmseqs easy-search $blstinfile $spdb $outfile tmp " 
                    . "--search-type 1 --format-mode 0 -e $eCutOff --threads $numProc -s 7.5 --max-seqs 500";
            }
            else {
                die "specified search method not supported. please use blastp, diamond, or mmseqs2\n";
            }
            
        my $blastOutput = `$blastCmd`;
        $t = localtime;
        print "$t complete blast $sp rc $blastOutput\n" if $verbose > 1;
        open( XX, ">>$OID_BLASTDIR/.$sp.Part.done" );  #create a done file
            print XX ".$sp.Part$prefix.done\n";
        close XX;

        #        my @blastResults = split /\n/, $blastOutput;
        #        foreach my $res(@blastResults){
        #            my $locus = (split /\t/, $res)[0];
        #            $blastRes{$locus} .= "$res\n";#$blastOutput;
        #        }
    }

    #    $count++;
    #    unlink $tmpFile;
    print "\n" if $verbose >= 1;

    #untie %blastRes;
    #untie %geneLen;

}

## rdblastres - reads in new blastres file/genlen file
sub rdblastrs {    #pass ref to %blastRes and %geneLen
    my $pblst      = shift;
    my $pgene      = shift;           # get 2 args
    my $blastResDB = $BLAST_RES_DB;
    my $geneLenDB  = $GENE_LEN_DB;
    my @stuff;
    my $count = 0;
    open FA, "$geneLenDB" or die "cannot open genelendb $geneLenDB\n";
    while (<FA>) {
        chomp;
        next if ( !$_ );
        @stuff = split;
        $$pgene{ $stuff[0] } = $stuff[1];
        $count++;
        print "gene $stuff[0] has genelen $$pgene{$stuff[0]}\n" if $count < 6;
    }
    print "loaded $count gene lengths\n";
    close FA;

    open FA, "$blastResDB" or die "cannot open blastresdb $blastResDB\n";
    my $locus = "";
    while (<FA>) {
        next if (/\s+\n/);
        if (/^>/) {
            if ($locus) {
                $count++;
                print "$locus\n$$pblst{$locus}" if ( $count < 3 );
            }
            chomp;
            $locus = substr( $_, 1 );
            $$pblst{$locus} = "" if ( !defined $$pblst{locus} );
            next;
        }
        $$pblst{$locus} .= $_;
    }    ## while fa
    close FA;
    print "loaded $count gene locus\n";
}

sub rdgeners {    #pass ref to blastres fa and %geneLen
    my $pblst      = shift;
    my $pgene      = shift;           # get 2 args
    my $blastResDB = $BLAST_RES_DB;
    my $geneLenDB  = $GENE_LEN_DB;
    my @stuff;
    my $count = 0;
    open FA, "$geneLenDB" or die "cannot open genelendb $geneLenDB\n";
    while (<FA>) {
        chomp;
        next if ( !$_ );
        @stuff = split;
        $$pgene{ $stuff[0] } = $stuff[1];
        $count++;
        print "gene $stuff[0] has genelen $$pgene{$stuff[0]}\n" if $count < 6;
    }
    print "loaded $count gene lengths\n";
    close FA;

    my $fa;
    open $fa, "$blastResDB" or die "cannot open blastresdb $blastResDB\n";
    $$pblst = $fa;

}

sub nextblst {    #get next gene - pass file id, array ref
    my $fa        = shift;
    my $pblst     = shift;
    my $plastgene = shift;
    my $lastgene  = $$plastgene;

    my @stuff;
    @$pblst = ();    #empty
    if ( !$lastgene ) {
        while (<$fa>) {
            next if (/\s+\n/);
            if (/^>/) {
                chomp;
                $lastgene = substr( $_, 1 );    #skip >
                print "first gene $lastgene\n";
                last;
            }
        }
    }
    return $lastgene if !$lastgene;
    while (<$fa>) {
        next if (/\s+\n/);
        if (/^>/) {
            chomp;
            $$plastgene = substr( $_, 1 );    # save next gene
            return $lastgene;
            last;
        }
        push @$pblst, $_;
    }
    $$plastgene = '';    #eof
    return $lastgene;
}

#
# MCL Clustering
# Edge weight is -log(e-value)*(alignment_ratio)^2
#

sub mclmcx() {
    my $log10 = log(10);
    my %geneLen;
    my $numProc    = $NCPU;
    my $mciFile    = "$OID_BLASTDIR/weights.mci";
    my $tabFile    = "$OID_BLASTDIR/weights.tab";
    my $blastResDB = $BLAST_RES_DB;
    my $gz;                            # gzip object
    my $thisFunc = ( caller(0) )[3];
    my $blfa     = 0;                  #fa for blastResDB
    my $lastgene = '';                 #keepsake for last fa read for blastres
    my $ct;

    $ct = localtime;
    print "$ct Generating edge weights for clustering ...\n";
    return if ( -f $mcxdone );         #don't redo
                                       #if ($DEBUG) {
     #    $gz = new IO::Compress::Gzip $weightFile or die "IO::Compress::Gzip failed :$GzipError\n";
     #}
     # Convert directly into undirected graph by using maximum of two edge weights for the same edge
     # and output info mcl matrix format.
    my $lsb = `ls -s $blastResDB`;
    $lsb =~ s/^\s//;
    my @bls   = split /\s/, $lsb;
    my $bllen = int( $bls[0] / 1000 );

    if ( $bllen < 1000 ) {    #gigabyte
        print "blastres len $bllen MB use mcxload stream-mirror\n";
        open( MFH,
"|mcxload -abc - --stream-mirror -re max -o $mciFile -write-tab $tabFile"
        );
    }
    else {
        print "big blastres len $bllen MB use mcxload ri max\n";

        open( MFH,
"|mcxload -abc - --write-binary -ri max -re max -o $mciFile -write-tab $tabFile"
        );
    }

    #    tie (%blastRes, "AnyDBM_File", $BLAST_RES_DB,O_RDONLY,0666)
    #      or die "Cannot open $BLAST_RES_DB: $!\n";
    #    tie (%geneLen, "AnyDBM_File", $GENE_LEN_DB,O_RDONLY,0666)
    #      or die "Cannot open $GENE_LEN_DB: $!\n";
    #    rdblastrs(\%blastRes,\%geneLen);
    rdgeners( \$blfa, \%geneLen );    # read in genlen db and open blastResDb
    while (1) {
        my @blastOutput = ();                                             #empty
        my $g           = nextblst( $blfa, \@blastOutput, \$lastgene );
        last if !$g;
        my %seen;
        my $seqLen = $geneLen{$g};
        if ( !$seqLen ) {
            print "no genelen entry for $g\n";
            next;
        }

        #        my @blastOutput = split( /\n/, $blastRes{$g});
        foreach (@blastOutput) {
            my ( $target, $alignLen, $eVal ) = (split)[ 1, 3, 10 ];

            # only use the first entry of each target
            next if defined $seen{$target};

            $seen{$target} = 1;
            my $targetSeqLen = $geneLen{$target};
            if ( !$targetSeqLen ) {
                print "no target genelen entry for $target\n";
                next;
            }
            my $alignRatio = $alignLen /
              ( ( $seqLen > $targetSeqLen ) ? $seqLen : $targetSeqLen );
            my $weight;
            if ( $eVal != 0 ) {
                $weight = -log($eVal) / $log10;
            }
            else {
                $weight = 500;
            }
            $weight *= $alignRatio * $alignRatio;

            #$gz->print("$g\t$target\t$weight\n") if $DEBUG;
            print MFH "$g\t$target\t$weight\n";
        }

    }

    #    untie %blastRes;
    #untie %geneLen;
    close MFH;

    #$gz->close() if $DEBUG;
    open( XX, ">$mcxdone" );    #create a done file
    close XX;
    $ct = localtime;
    print "$ct mcxload done\n";
}

sub mclCluster() {
    my $mcl_exe = "mcl";
    my %lengthOfGene;
    my $log10 = log(10);
    my %blastRes;
    my %geneLen;
    my $numProc    = $NCPU;
    my $mciFile    = "$OID_BLASTDIR/weights.mci";
    my $tabFile    = "$OID_BLASTDIR/weights.tab";
    my $weightFile = "$OID_BLASTDIR/weights.dat.gz";    # for debugging only
    my $blastResDB = $BLAST_RES_DB;
    my $gz;                                             # gzip object
    my $thisFunc = ( caller(0) )[3];
    my $blfa     = 0;                  #fa for blastResDB
    my $lastgene = '';                 #keepsake for last fa read for blastres
    my $ct;

    $ct = localtime;
    print "$ct Generating edge weights for clustering ...\n";
    if ( !-f $mcxdone ) {

#if ($DEBUG) {
#    $gz = new IO::Compress::Gzip $weightFile or die "IO::Compress::Gzip failed :$GzipError\n";
#}
# Convert directly into undirected graph by using maximum of two edge weights for the same edge
# and output info mcl matrix format.
        my $lsb = `ls -s $blastResDB`;
        $lsb =~ s/^\s//;
        my @bls   = split /\s/, $lsb;
        my $bllen = int( $bls[0] / 1000 );
        if ( $bllen < 1000 ) {    #gigabyte
            print "blastres len $bllen MB use mcxload stream-mirror\n";
            open( MFH,
"|mcxload -abc - --stream-mirror -re max -o $mciFile -write-tab $tabFile"
            );
        }
        else {
            print "big blastres len $bllen MB use mcxload ri max\n";

            open( MFH,
"|mcxload -abc - --write-binary -ri max -re max -o $mciFile -write-tab $tabFile"
            );
        }

        #    tie (%blastRes, "AnyDBM_File", $BLAST_RES_DB,O_RDONLY,0666)
        #      or die "Cannot open $BLAST_RES_DB: $!\n";
        #    tie (%geneLen, "AnyDBM_File", $GENE_LEN_DB,O_RDONLY,0666)
        #      or die "Cannot open $GENE_LEN_DB: $!\n";
        #    rdblastrs(\%blastRes,\%geneLen);
        rdgeners( \$blfa, \%geneLen );   # read in genlen db and open blastResDb

        #    foreach my $g ( keys %blastRes ) {
        while (1) {
            my @blastOutput = ();        #empty
            my $g           = nextblst( $blfa, \@blastOutput, \$lastgene );
            last if !$g;
            my %seen;
            my $seqLen = $geneLen{$g};
            if ( !$seqLen ) {
                print "no genelen entry for $g\n";
                next;
            }

            #        my @blastOutput = split( /\n/, $blastRes{$g});
            foreach (@blastOutput) {
                my ( $target, $alignLen, $eVal ) = (split)[ 1, 3, 10 ];

                # only use the first entry of each target
                next if defined $seen{$target};

                $seen{$target} = 1;
                my $targetSeqLen = $geneLen{$target};
                if ( !$targetSeqLen ) {
                    print "no target genelen entry for $target\n";
                    next;
                }
                my $alignRatio = $alignLen /
                  ( ( $seqLen > $targetSeqLen ) ? $seqLen : $targetSeqLen );
                my $weight;
                if ( $eVal != 0 ) {
                    $weight = -log($eVal) / $log10;
                }
                else {
                    $weight = 500;
                }
                $weight *= $alignRatio * $alignRatio;

                #$gz->print("$g\t$target\t$weight\n") if $DEBUG;
                print MFH "$g\t$target\t$weight\n";
            }

        }

        #    untie %blastRes;
        #untie %geneLen;
        close MFH;
    }
    else {
        print "mcxload already done\n";
    }    # if not mcxdone
         #$gz->close() if $DEBUG;
    open( XX, ">$mcxdone" );    #create a done file
    close XX;

    # run mcl
    $ct = localtime;
    print "$ct Running mcl ...\n";

    #my @mclArgs = ($mclInput, "--abc", "-I", "1.4", "-o", $mclOutput);
    my @mclArgs = (
        $mciFile,   "-force-connected", "y",  "-scheme", "7", "-I", "1.4",
        "-use-tab", $tabFile,           "-o", $clustersFile
    );

    #    push(@mclArgs, "-te", $numProc) if $numProc > 1;
    my $status = system( $mcl_exe, @mclArgs );
    die "$thisFunc: $mcl_exe exit status = $?" unless $status == 0;
    $ct = localtime;
    print "$ct mcl done.\n";
    open( XX, ">$clusterdone" );    #create a done file
    close XX;
}

sub make1fam {    # makes a single family
    my $verbose   = $VERBOSE;
    my $familyNum = shift;
    my $clust     = shift;
    my $pfasta    = shift;
    my $familyDir;

#    my $FASTACMD = defined($BLAST_HOME) ? "$BLAST_HOME/bin/fastacmd" : "fastacmd";
    my $FASTACMD = "blastdbcmd";

    my $seq;
    my @out;
    my @genes = split /\s/, $clust;
    if ( @genes == 1 ) {
        return (1);
    }
    else {
        print "Writing family $familyNum (" . scalar(@genes) . " sequences)\n"
          if $verbose == 2;
        $familyDir = "$OID_DATADIR/$familyNum";
        if ( !-d $familyDir ) {
            mkdir $familyDir, 0755 or die "Cannot create $familyDir: $!";
        }
        open FFH, ">$familyDir/$unalignedFamily"
          or die "Cannot open $familyDir/$unalignedFamily for writing: $!";
        my %spsInFam;
        foreach my $gene (@genes) {
            if ( !exists $$pfasta{$gene} ) {
                print "$gene not in fasta table\n";
                next;
            }
            print FFH ">$gene\n";
            print FFH $$pfasta{$gene};    #has imbeded carrage returns
        }

    }
    close FFH;
    open( XX, ">>$OID_DATADIR/.fam.done" );    #add to done file
        print XX ".$familyNum.fam.done\n";
    close XX;
    my $ngene = scalar(@genes);
    return ($ngene);
}    #make1fam

sub makeFamily() {
    my $verbose    = $VERBOSE;
    my $familyNum  = 1;
    my $singletNum = 0;
    my $familyDir;
    my $singletFile = "$OID_DATADIR/singlets";
    my $status;

#    my $FASTACMD = defined($BLAST_HOME) ? "$BLAST_HOME/bin/fastacmd" : "fastacmd";
    my $FASTACMD = "blastdbcmd";

    # Clustering
    setglsttm();
    mclCluster if !-f $clusterdone;
    if ( -f $familydone ) {
        print "families already done\n";
        return;
    }

    # Output gene families
    open MFH, $clustersFile or die "Cannot open $clustersFile: $!\n";

    open SFH, ">>$singletFile"
      or die "Create open $singletFile for writing: $!\n";

    while (<MFH>) {
        chomp;
        if ( system( "grep -q '\\.$familyNum\\.fam\\.done' $OID_DATADIR/.fam.done") < 1 ) {
            $familyNum++;
            next;
        }
        my $seq;
        my @out;
        my @genes = split;
        if ( @genes == 1 ) {
            print SFH "$genes[0]\n";
            $singletNum++;
        }
        else {
            print "Writing family $familyNum ("
              . scalar(@genes)
              . " sequences)\n"
              if $verbose == 2;
            $familyDir = "$OID_DATADIR/$familyNum";
            if ( !-d $familyDir ) {
                mkdir $familyDir, 0755 or die "Cannot create $familyDir: $!";
            }
            open FFH, ">$familyDir/$unalignedFamily"
              or die "Cannot open $familyDir/$unalignedFamily for writing: $!";
            my %spsInFam;
            foreach my $gene (@genes) {
                my $sp = ( parseOIDgeneName($gene) )[0];
                $spsInFam{$sp} .= "$gene\n";
            }
            foreach my $sps ( keys %spsInFam ) {
                open TMPIDS, ">tempIDs";

                #$spsInFam{$sps} =~ s/\s/\n/g;
                print TMPIDS $spsInFam{$sps};
                close TMPIDS;
                @out =
                  `$FASTACMD -db $sps -entry_batch tempIDs -outfmt '>%a %s'`;
                for ( my $i = 0 ; $i <= $#out ; $i++ ) { $out[$i] =~ s/\s/\n/; }
                print FFH @out;
            }

            #            foreach my $gene (@genes) {
            #                my $sp = (parseOIDgeneName($gene))[0];
            #                if ( grep /^$sp$/, ( @INGROUP, @OUTGROUP ) ) {
            #                    @out = `$FASTACMD -d $sp -s $gene`;
            #                    @out = `$FASTACMD -db $sp -entry $gene`;
            #                }
            #                if (@out == 0) {
            #                    warn "Failed to look up $gene in \"$sp\": $?";
            #                    next;
            #                }
            #                shift @out;
            #                my $seq = join "", @out;
            #                print FFH ">$gene\n$seq\n";
            #            }
            close FFH;
            open( XX, ">$OID_DATADIR/.fam.done" )
              ;    #create a done file
                print XX ".$familyNum.fam.done\n";
            close XX;
            $status = system("$OID_HOME/bin/fix_familyfastas.sh $familyDir/$unalignedFamily"); #fix SEQUENCE>TaxID errors
            chkuptm();
            $familyNum++;
        }
    }

    close MFH;
    open( XX, ">$familydone" );    #create a done file
    close XX;

    print $familyNum - 1 . " families\n";
    print "$singletNum singlets\n";
}

#
# Create alignment for each family (with elision).  Need algn_elide.sh.
# Optional arg: family num pattern
#
sub alignFamily {

    my $dirRE = shift;
    if ( !defined($dirRE) ) {
        print "aligned family called withot argument\n";
        return;
    }

    #    $dirRE = '.*' if !defined($dirRE);

    my $dir = $dirRE;
    if ( $dirRE =~ /\^(\d+)\$/ ) {
        $dir = $1;
    }

    #    print "$dir, $dirRE\n";
    my $verbose = $VERBOSE;

    # Defaults
    my $thisFunc = ( caller(0) )[3];

    # Change to data dir
    my $savDir = getcwd;
    chdir $OID_DATADIR;

    my $count = 0;
    print "Generating alignments ...\n" if $verbose;

    # Go over each family
    #    foreach my $dir (<[1-9]*>) {
    #        next if $dir !~ /$dirRE/;
    $count++;
    chdir $dir
      or die "$thisFunc: failed to change to family directory $dir: $!\n";
    if ( -s $alignedFamily ) {
        print "Skipping family $dir - alignment file exists\n"
          if $verbose == 2;

        #            chdir $OID_DATADIR;
        chdir $savDir;

        #            next;
        return;
    }
    else {
        print "Creating alignment for family $dir\n" if $verbose == 2;
    }

    my $status;

    # Aligning
    if ( $verbose == 2 ) {
        $status = system("align_family.sh $unalignedFamily $alignedFamily");
    }
    else {
        if ( $verbose == 1 ) {
            $| = 1;    # flush STDOUT
            print "\b" x 7;
            printf "%6d", $count;
            $| = 0;
        }
        $status = system(
            "align_family.sh $unalignedFamily $alignedFamily >/dev/null 2>&1");
    }
    warn "$thisFunc: alignment error - exit status = $?"
      unless $status == 0;

    #        chdir $OID_DATADIR;
    #    }

    print "\n" if $verbose == 1;
    chdir $savDir;
}

#
# Create guide tree for each family.
# Optional arg: family num pattern
#
sub makeTree {

    my $dirRE = shift;
    if ( !defined($dirRE) ) {
        print "aligned family called withot argument\n";
        return;
    }

    #    $dirRE = '.*' if !defined($dirRE);
    my $dir = $dirRE;
    if ( $dirRE =~ /\^(\d+)\$/ ) {
        $dir = $1;
    }
    my $verbose = $VERBOSE;

    # Defaults
    my $prefix = "oid";

    #my $paup_exe  = "paup";
    my $tnt_exe      = "tnt";
    my $nrat         = 10;       # number of ratchets
    my $niter        = 100;      # number of iterations per ratchet
    my $timeLimit    = 10800;    # Time limit for a heuristic search
    my $ratTimeLimit = 60;  # Time limit for a search within a ratchet iteration

    my $thisFunc = ( caller(0) )[3];
    my $oldFH;
    my %seq;
    my @outgroup;           # Outgroup taxa for current family
    my $curcut  = 0;
    my $maxcut  = 0;
    my $curproc = 0;
    my $ntaxa   = 0;

    # Change to data dir
    my $savDir = getcwd;
    chdir $OID_DATADIR;

    # Go over each family
    my $count = 0;
    print "Making trees ...\n" if $verbose;

    #    foreach my $dir (<[1-9]*>) {
    #        next if $dir !~ /$dirRE/;
    for ( my $x = 0 ; $x < 1 ; $x++ ) {
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
            print "$dir Tree file already exists ... skipping\n"
              if $verbose;

            #            chdir $OID_DATADIR;
            #            next;
            chdir $savDir;
            return;
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
        if ( $ntaxa < 4 ) {
            print
"$dir Too few taxa. TNT cannot create tree with less than 4 taxa\n"
              if $verbose;

            #            chdir $OID_DATADIR;
            #            next;
            chdir $savDir;
            return;
        }

        #        if (-s $prefix.".nex") goto RUNTNT
        open( NEX, ">", $prefix . ".nex" )
          or die "$thisFunc: Failed to open file for writing: $!\n";
        $oldFH = select(NEX);
        print "xread\n";
        print "$length $ntaxa\n";
        print "&[proteins nogaps]\n";
        foreach ( sort keys %seq ) {
            if ( length($_) < 32 ) {

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
        #
        #     we find proc to write by cutoff
        #
        for ( my $i = 0 ; $i <= $#cutoff ; $i++ ) {
            if ( $ntaxa <= $cutoff[$i] ) {
                $curcut = $cutoff[$i];
                $maxcut = $cutmax[$i];
                last;
            }
        }
        if ( !$curcut ) {
            $curcut = $cutoff[$#cutoff];
            $maxcut = $cutmax[$#cutoff];
            print
"$ntaxa greater that largest cutoff $curcut max $maxcut - family probably too large\n";
        }
        if ( $maxcut < 2 ) {    #only 1 proc
            $curproc = 1;
            open( PROC, ">$prefix" . ".proc" )
              or die "$thisFunc: Failed to open file for writing: $!\n";
            $oldFH = select(PROC);

            # TNT commands from %procptr
            my $ar = $procptr{$curcut}->[$curproc];
            foreach my $line (@$ar) {
                print "$line\n";
            }
            close(PROC);

            select($oldFH);

            # Execute tnt

            ###RUNTNT:
            my $status = system("tnt p $prefix.proc 0</dev/null 1>&0 2>&0");
            die "$thisFunc: TNT error: $?" unless $status == 0;

            #            print "would have run tnt p $prefix.proc\n";
        }
        else {    #here maxproc > 1 so we must do a multiprocess
            for ( $curproc = 1 ; $curproc <= $maxcut ; $curproc++ ) {
                open( PROC, ">$prefix$curproc" . ".proc" )
                  or die "$thisFunc: Failed to open file for writing: $!\n";
                $oldFH = select(PROC);

                #     TNT commands from %procptr
                my $ar = $procptr{$curcut}->[$curproc];
                foreach my $line (@$ar) {
                    print "$line\n";
                }
                close(PROC);
                select($oldFH);
            }
            my $status = system("runtnta.pl $maxcut ");

            #           print "would have run runtnt.py $prefix $maxcut\n";
        }
        unlink <$prefix.tre.*>;
        unlink "$prefix.tmp" if -f "$prefix.tmp";

        #            chdir $OID_DATADIR;
    }

    print "\n" if $verbose == 1;
    chdir $savDir;

}

sub restartme {

    #    my $PARAMS = "-l nodes=1:ppn=1 ";
    #    my $MY_SCRIPT="qsub $PARAMS $OID_WRAPPER $OID_USER_DIR/pipe.pbs";
    my $MY_SCRIPT = "$OID_HOME/bin/pipe.pbs";
    if ( $ENV{"MY_SHELL"} ) {
        $MY_SCRIPT = $ENV{"MY_SHELL"};
    }

    #    my $myjob = $ENV{"PBS_JOBID"};
    #    my $submit = q/qsub  /."$PARAMS $MY_SCRIPT";
    my $submit = "";
    if ( $HPC =~ /S/ ) {
        $MY_SCRIPT = "-o toplog/%J.out $OID_WRAPPER $MY_SCRIPT";
        $submit    = q/sbatch / . "$MY_SCRIPT";
    }
    else {
        $submit = q/qsub / . "$MY_SCRIPT";
    }
    print "$submit\n";
    $! = 1;    #flush
    my $pid = open( QS, "$submit |" );
    my $qs  = <QS>;
    if ($qs) {
        chomp $qs;
        print "start qsid $qs for $MY_SCRIPT\n";
        print "I am restarting myself\n";
    }
    print "I mush cancel my job pid $myjob\n";
    $| = 1;    #flush
    if ( $HPC =~ /S/ ) {
        my $rc = system("scancel $myjob");
    }
    else {
        my $rc = system("qdel $myjob");
    }
    sleep 5;
    print "I should not get here\n";
    $| = 1;     #flush
    exit(0);    ## exit
}    #restarts me

sub setlimit {    #format string nnGB,walltime=mm:00:00"
    my $fam   = shift;                #must have a arg
    my $fammx = $limkey[$#limkey];    #default falls through loop
    foreach my $val (@limkey) {
        if ( $fam < $val ) {
            $fammx = $val;
            last;
        }
    }
    my $set1   = $limval{$fammx}[0];
    my $set2   = $limval{$fammx}[1];
    my $setstr = "mem=$set1" . "GB,walltime=" . "$set2" . ":00:00";
    $NCPU = $cutmax[$#cutoff];
    for ( my $i = 0 ; $i <= $#cutoff ; $i++ ) {
        if ( $fam <= $cutoff[$i] ) {
            $NCPU = $cutmax[$i];
            last;
        }
    }
    print "setlimit set $setstr  ncpu $NCPU\n";
    return $setstr;
}

sub strtgrp {    # starts group qsub
    my $ngroup    = shift;
    my $grouplo   = shift;
    my $grouphi   = shift;
    my $GP_SCRIPT = "$OID_HOME/bin/run_tgroup" . $HPC . ".sh";
    my $grpcpu    = 20;
    if ( $ngroup < 40 ) {
        $grpcpu = int( $ngroup / 2 );
        $grpcpu = 2 if $grpcpu < 2;
    }
    my $gmem = 4 * $grpcpu;
    print "group cpu $grpcpu mem $gmem", "GB\n";

#        my $PARAMS = "-l nodes=1:ppn=$grpcpu -l mem=$gmem"."GB,walltime=12:00:00";
#        my $submit = q/qsub  /."$PARAMS $OID_WRAPPER $GP_SCRIPT -v arg1=$grouplo,arg2=$grouphi,arg3=$grpcpu".' 2>/dev/null';
    my $nt   = 1;
    my $wall = "12:00:00";
    my $gb   = "$gmem" . "GB";
    my @args = ();
    push @args, $grouplo;
    push @args, $grouphi;
    push @args, $grpcpu;

    my $submit = fmthpc( $gb, $wall, $grpcpu, $nt, $GP_SCRIPT, @args );
    $submit .= " 2>/dev/null";
    print "$submit\n";
    my $pid = open( QS, "$submit |" );
    my $qs  = <QS>;
    chomp $qs if $qs;
    return fixqs($qs);
}

sub strbmon {    # starts Bigmon
    my $GP_SCRIPT = "$OID_HOME/bin/run_bigmon" . $HPC . ".sh";
    #        my $PARAMS = "-l nodes=1:ppn=1 -l mem=12"."GB,walltime=24:00:00";
    #        my $submit = q/qsub  /. "$PARAMS $GP_SCRIPT ".' 2>/dev/null';
    my $nt   = 1;
    my $wall = "24:00:00";
    my $gb   = "12GB";
    my $mcpu = 1;
    my @args = ();
    if ( $TREEPROGRAM =~ /OBLONG/ ) {
        $GP_SCRIPT = "-o $OID_USER_DIR/log/job/bigmon.log $OID_HOME/bin/Bigmon2.sh";
        @args = ($MAXQS);
    }
    my $submit = fmthpc( $gb, $wall, $mcpu, $nt, $GP_SCRIPT, @args );
    print "$submit\n";
    my $pid = open( QS, "$submit |" );
    my $qs  = <QS>;
    chomp $qs if $qs;
    return fixqs($qs);
}

sub startqs {    #starts a qsub

    #                startqs(\@qsid,$cmpqs,$famx);
    ( my $rqs, my $cmpqs, my $fam ) = @_;
    my $walmem     = setlimit( $sztaxa{$fam} );
    my $JOB_SCRIPT = "$OID_USER_DIR/run_oid_job.sh";
    my $PARAMS     = "-l nodes=1:ppn=$NCPU -l $walmem";

    #    my $PARAMS = "-l nodes=1:ppn=$NCPU -l mem=8GB,walltime=12:00:00";
    my $submit =
        q/qsub  /
      . "$PARAMS $JOB_SCRIPT"
      . q/ -v arg1="-at",arg2="^/ . "$fam"
      . '$" 2>/dev/null';
    print "$submit\n";

#    my $pid = open(QS,"qsub -l $PARAMS $JOB_SCRIPT -v arg1=/"$at/",arg2=/"^$fam.$dol/" 2>/dev/null |");
    my $pid = open( QS, "$submit |" );

    my $qs = <QS>;
    if ($qs) {
        chomp $qs;
        print "start qsid $qs for family $fam\n";
        $$rqs{$qs} = $fam;    #put $qs in its place
        return 1;
    }
    else {
        print "qsub somehow fails\n";
        ## bug in code - will never end
        return 0;
    }
}    ## startqs

sub activqs {    #test if qsid still running

#modified to only check for group and bigmon - we dont qsid for big taxa anymore
    my $rqs      = shift;
    my $qsgrp    = shift;
    my $qsbigmon = shift;    #new task
    my %active   = ();
    my @stuff    = ();
    my $qs;
    my $fam;
    my %actfam = ();

    #    my $myjob = $ENV{"PBS_JOBID"};
    my $fndme = 0;
    my $flag  = 1;
    $fndme = 1 if ( !$myjob );    #non pbs system
    if ( !$myjob ) {              #not running in qsid
        $fndme = 1;
        $myjob = -1;              #to avoid error msg
    }

    #    print "testqs called @$rqs\n";
    $| = 1;                       # flush STDOUT
    getqids( \%active );
    if ($myjob) {
        $fndme = 1 if exists $active{"$myjob"};
    }

#    print "enter testqs as $MYUSER\n";
#    while ($flag) {
#        $flag = 0;
#    my $pid = open(QS,"qstat -e -u $MYUSER |");
#    my $jnk = <QS>;
#    $jnk = <QS>;
#    print "testqs $jnk\n";
#    while (<QS>){
#        chomp;
#        @stuff = split;
#        next if ($stuff[0] =~ /\D/);
#        print "qs $stuff[0], $stuff[-2]\n";
#        now that we are not looking for individual tnt task - don't check if not started yet
#        if ($stuff[-2] =~ /R|Q/){
#            if ($stuff[-2] =~ /Q/){
#                print "$stuff[0] still on Q wait for start\n";
#                $flag = 1;
#                sleep 30;
#                last;
#            }
#            $qs = "$stuff[0]";

    #            $fndme = 1 if ($qs == $myjob);
    #            next if ($qs == $myjob);
    #            $active{"$qs"} = $stuff[-1];  #put found pids in hash
    #            print "$stuff[0] is active\n";
    #        }
    #    }
    #} ## while flag
    if ( !$fndme ) {
        print "did not find me job $myjob - something wrong - try later\n";
        return ( $fndme, %actfam );
    }
    foreach $qs ( keys %active ) {

        #        my $pid = open(GR, "grep -c '>' FAMILY |") ;
        #        $sztaxa{$dir} = <GR>;
        my $log = "log/job/*" . $qs;

        #        print $log
        #        don't check for running jobs anymore
        #        my $pid = open(GR, "grep 'Tree search' $log 2>/dev/null| ");
        #        $fam = 0;
        #        $_ = <GR>;
        #        next if (! $_);
        #        chomp;
        #        print "found $_\n";
        #        if (/\D+(\d+)/){
        #            $fam = $1;
        #        }
        #        print "active $qs found fam $fam wall $active{$qs}\n";
        #        $actfam{$fam} = $qs if ($fam>0);  #reverse hash
        #        $$rqs{$qs} = $fam  if ($fam>0);
        if ( $qs == $qsgrp ) {
            print "qsid $qsgrp active\n";
            $$rqs{$qs} = 1;
        }
        if ( $qs == $qsbigmon ) {
            print "qsid $qsbigmon active\n";
            $$rqs{$qs} = 1;
        }
    }
    return ( $fndme, %actfam );
}

sub testqs {    #test if qsid still running
    my $rqs    = shift;
    my %active = ();
    my @stuff  = ();
    my $qs;
    my $fndme = 0;
    $fndme = 1 if ( !$myjob );    #non pbs system

    #    print "testqs called @$rqs\n";
    $| = 1;                       # flush STDOUT
    getqids( \%active );

    #    my @kq = keys %active;
    #    print "keys active ",join ",",@kq,"\n";
    #    print "enter testqs as $MYUSER\n";
    if ($myjob) {
        $fndme = 1 if exists $active{"$myjob"};
    }
    if ( !$fndme ) {
        print "did not find me job $myjob - something wrong - try later\n";
        return ( 0, 0 );
    }
    for $qs ( keys %$rqs ) {

        #        print "test active $qs\n";
        next if ( $qs =~ /D/ );

        #        print "$qs is still active\n" if exists $active{"$qs"};
        next if exists $active{$qs};

        #        foreach my $x (%active){
        #            print "$x not number\n " if ($x =~/D/);
        #            print "$x and $active{$x}\n";
        #            print "$qs is a keys\n " if ($qs == $x);
        #        }
        print "really keys active ", join ":", keys %active, "\n";

        #        my @blinfo = @{$$rqs{$qs}};

        my $blinfo = $$rqs{$qs};

      #        print "$qs  blast $blinfo not active\n" if ! exists $active{$qs};
        delete $$rqs{$qs};

        #        return ($qs,@blinfo);
        return ( $qs, $blinfo );
    }
    return ( 0, 0 );    #everthing still running or qued
}    #end testqs

sub schedatqs {
    my $restfl = shift;    #to allow a restart with running qsubs
    our $familyDB = "$OID_USER_DIR/familyDB.dbm";
    my $dol    = '$';
    my $dupcnt = 0;

    #my %sztaxa = ();
    my $grouplo = 0;       #lowest dir to run in group
    my $grouphi = 0
      ; #highest dir to run in group -- note grouplo has larger family size (lower group, bigger family
    my $ngroup = 0
      ; #count of directories in group (may be missing dir between grouplo and grouphi
    my $nbigmon  = 0;   #count of directories used by bigmon - just for interest
    my $tgroupsz = 0;   #cum size of all family in group - just for interest

    my @thisdis = ( 0, 0, 0, 0, 0 );
    my $tstrt   = localtime;
    print "start time $tstrt\n";
    print "in scedatqs arg $restfl\n" if $restfl;
    print "my user is $MYUSER\n";
    setglsttm();

    #    my $t = time();
    #    $tstrt = $t;
    #    if (-s 'starttime'){
    #        open FI,"<starttime";
    #        my $x = <FI>;
    #        close FI;
    #        if ($x){
    #            chomp($x);
    #            my $dif = $t - $x;
    #            my $xx = localtime($x);
    #            print "$xx in starttime used to time us, dif $dif \n";
    #            $tstrt = $x;
    #        }
    #    }
    my $savdir = getcwd;
    chdir $OID_DATADIR;
    unlink "$SCHEDONE" if ( -f "$SCHEDONE" );    #delete if exists
    my $badtie = 0;
    tie( %sztaxa, "DB_File", $familyDB, O_RDONLY, 0666 )
      or do { print "could not create tie $familyDB\n"; $badtie = 1; };
    if ($badtie) {
        foreach my $dir (<[1-9]*>) {
            chdir $dir;
            my $pid = open( GR, "grep -c '>' FAMILY |" );
            $sztaxa{$dir} = <GR>;
            chdir $OID_DATADIR;
        }
    }
    my @kefnd = keys %sztaxa;
    my $fam   = scalar @kefnd;
    print "there are $fam families\n";

    #    my @taxsz = sort { $sztaxa{$a} <=> $sztaxa{$b} } @kefnd;
    my @taxsz = sort { $a <=> $b } @kefnd;
    for $fam (@taxsz) {
        my $thissz = $sztaxa{$fam};
        if ( $thissz > $TNTA and $thissz < $TNTB ) {
            $grouplo = $fam if !$grouplo;    #first family in range
            $grouphi = $fam;
            $ngroup++;
            $tgroupsz += $thissz;
        }
        $nbigmon += 1
          if $thissz > $TNTB;    # this is # of direcories sched by Bigmon

        if ( $thissz < 50 ) {
            $thisdis[0]++;
        }
        elsif ( $thissz < $TNTA ) {
            $thisdis[1]++;
        }
        elsif ( $thissz < $TNTB ) {
            $thisdis[2]++;

            #            print "$fam has size $thissz\n";
        }
        elsif ( $thissz < $TNTC ) {
            $thisdis[3]++;
            print "$fam has size $thissz\n";
        }
        else {
            $thisdis[4]++;
            print "$fam has size $thissz\n";
        }
    }
    print "counts @thisdis\n";
    print "<50,$TNTA,$TNTB,$TNTC,big\n";
    chdir $savdir;
    my @famtbl   = reverse @taxsz;
    my $qsct     = 0;
    my $totqs    = 0;
    my $lastqs   = 0;
    my $lastbig  = 0;
    my $qsgrp    = 0;
    my $qsbigmon = 0;
    my %qsid     = ();
    my $oldqs;
    my $qsfam      = 0;
    my $schedrst   = 0;
    my %actfam     = ();
    my $bigest     = 999999;
    my $JOB_SCRIPT = "$OID_USER_DIR/run_oid_job.sh";

    print "tree program $TREEPROGRAM \n";
    if ( $TREEPROGRAM =~ /OBLONG/ ) {
        if ( -s "$SCHEDLOG" ) {
            open( SQ, "$SCHEDLOG" );
            while (<SQ>) {
                chomp;
                next if (/^\s+/);
                my $pid = $_;
                if (/^bigmon (\d+)/) {
                    $qsbigmon = $1;
                    next;
                }
            }
            close SQ;
        }
        else {
            unlink "$BMDONE";    ## first time sb removed
        }
        my $fndme = 0;
        until ($fndme) {         ## sometimes qstat doesn't work - so wont find me
            ( $fndme, %actfam ) = activqs( \%qsid, $qsgrp, $qsbigmon );
            sleep 15 if ( !$fndme );
        }
        $qsgrp = 99999999;
        $qsbigmon = 0 if ( $qsbigmon and !exists $qsid{$qsbigmon} );
        $qsct     = keys %qsid;    #init with actual active tasks
        print "were running $qsct tasks on restart\n";
   
        if ($ngroup) {
            $tgroupsz = $tgroupsz / $ngroup;
            print
    "there are $ngroup families between $grouplo and $grouphi avg family: $tgroupsz\n";
        }
        if ($nbigmon) {
            print "there are $nbigmon families with ntaxa > $TNTB\n;";
        }
        if ($qsbigmon) {
            print "$qsbigmon is Bigmon pid - still running\n";
        }
        elsif ( -f "$BMDONE" ) {
            print "Bigmon done\n";
            $qsbigmon = 9999999;
        }
        else {
            my $qs = strbmon();
            if ($qs) {
                print "start qsid $qs for  Bigmon\n";
                $qsid{$qs} = 1;
                $qsbigmon = $qs;
                $qsct++;
                $totqs++;
            }
            else {
                print "qsub failed\n";
            }
        }
        open( SQ, ">$SCHEDLOG" );
        print SQ "group $qsgrp\n";
        print SQ "bigmon $qsbigmon\n";
        close SQ;
        sleep 15;           #give a little time for  things to settle
        %qsid   = ();       #empty
        $fndme  = 0;
        %actfam = ();       #empty
        until ($fndme) {    ## sometimes qstat doesn't work - so wont find me
            ( $fndme, %actfam ) = activqs( \%qsid, $qsgrp, $qsbigmon );
        }
    
        $qsct = keys %qsid;    #init with actual active tasks
        print "$qsct qsid \n";
        my $tmnw = localtime;
        print "$tmnw\n";
        $lastbig = 1;
    
        my $flag = 0;
        $lastqs = 1;           # must keep 1 because bigmon does all qsubs
        while ( $qsct > 0 ) {
            ( $oldqs, $qsfam ) = testqs( \%qsid );
            $flag = 0;
            if ($oldqs) {
                print "complete fam $qsfam qcnt $qsct\n";
                $qsct--;
                $flag = 1;
                if ( $oldqs == $qsbigmon and !-f "$BMDONE" ) {
                    my $qs = strbmon();
                    if ($qs) {
                        $qsid{$qs} = 1;
                        $qsbigmon = $qs;
                        $qsct++;
                        print "start qsid $qs for  Bigmon\n";
                        open( SQ, ">$SCHEDLOG" );
                        print SQ "group $qsgrp\n";
                        print SQ "bigmon $qsbigmon\n";
                        close SQ;
                    }
                    else {
                        print "could not start Bigmon\n";
                    }
                    next;
                }
            }
            next if ($flag);
    
            ## later put in code to restart outselves if we are too old
            print "sleeping\n";
            sleep 300;
            chkuptm();
    
        }
    
        open( SQ, ">$SCHEDLOG" );
        print SQ "group $qsgrp\n";
        print SQ "bigmon $qsbigmon\n";
        close SQ;
        print "all $totqs qsubs ended\n";
        open( FI, ">$SCHEDONE" );
        print FI "$totqs\n";
        close FI;

    }
    else { #TNT
        print "tnt\n";
        if ( -s "$SCHEDLOG" ) {
            open( SQ, "$SCHEDLOG" );
            while (<SQ>) {
                chomp;
                next if (/^\s+/);
                my $pid = $_;
                if (/^group (\d+)/) {
                    $qsgrp = $1;
    
                    #                $qsid{$qsgrp} = 1;
                    #                $qsct++;
                    #                $totqs++;
                    next;
                }
                if (/^bigmon (\d+)/) {
                    $qsbigmon = $1;
                    next;
                }
    
                #            $qsid{$pid} = 1;
                #            $qsct++;
                #            $totqs++;
            }
            close SQ;
    
            #        for (my $i=0;$i<$qsct;$i++){
            #            ($oldqs,$qsfam) = testqs(\%qsid);
            #            if ($oldqs) {
            #                $qsct--;
            #                $qsgrp = 0 if ($oldqs == $qsgrp);
    ##            }
            #        }
    
        }
        else {
            unlink "$PLDONE";    ## first time sb removed
            unlink "$BMDONE";    ## first time sb removed
        }
        my $fndme = 0;
        until ($fndme) {         ## sometimes qstat doesn't work - so wont find me
            ( $fndme, %actfam ) = activqs( \%qsid, $qsgrp, $qsbigmon );
            sleep 15 if ( !$fndme );
        }
    
        #        $bigest = 999999;
        #        if (@actfam){
        #            @actfam = sort {$a <=> $b} @actfam;
        #            $bigest = $actfam[$#actfam];
        #        }
        $qsgrp    = 0 if ( $qsgrp    and !exists $qsid{$qsgrp} );
        $qsbigmon = 0 if ( $qsbigmon and !exists $qsid{$qsbigmon} );
        $qsct     = keys %qsid;    #init with actual active tasks
        print "were running $qsct tasks on restart\n";
    
        if ($ngroup) {
            $tgroupsz = $tgroupsz / $ngroup;
            print
    "there are $ngroup families between $grouplo and $grouphi avg family: $tgroupsz\n";
            if ($qsgrp) {
                print "$qsgrp is group pid - still running\n";
            }
            elsif ( -f "$PLDONE" ) {
                print "pool task complete\n";
                $qsgrp = 99999999;    #keeps things simpler
    
            }
            else {
                my $qs = strtgrp( $ngroup, $grouplo, $grouphi );
                if ($qs) {
                    print "start qsid $qs for  group\n";
                    $qsid{$qs} = 1;
                    $qsgrp = $qs;
                    $qsct++;
                    $totqs++;
                }
                else {
                    print "qsub failed\n";
                }
            }
        }
        if ($nbigmon) {
            print "there are $nbigmon families with ntaxa > $TNTB\n;";
            if ($qsbigmon) {
                print "$qsbigmon is Bigmon pid - still running\n";
            }
            elsif ( -f "$BMDONE" ) {
                print "Bigmon done\n";
                $qsbigmon = 9999999;
            }
            else {
                my $qs = strbmon();
                if ($qs) {
                    print "start qsid $qs for  Bigmon\n";
                    $qsid{$qs} = 1;
                    $qsbigmon = $qs;
                    $qsct++;
                    $totqs++;
                }
                else {
                    print "qsub failed\n";
                }
            }
        }
        open( SQ, ">$SCHEDLOG" );
        print SQ "group $qsgrp\n";
        print SQ "bigmon $qsbigmon\n";
    
        #    if (-s "schedrst") { #we restarted ourselves
        #       open (SQ,"schedrst");
        #        $_ = <SQ>;
        #        chomp;
        #        if (/\d+/){
        #            $lastbig = $_;
        #        }
        #        print "saw schedrst  at lastbig $lastbig = could do better job\n";
        #        $schedrst = 1 if $lastbig;
        #    }
    
        #    for $fam (@ntax)  #smaller fam is larger Family size
        $lastqs = 1;    #simple way to skip next loop without deleting a lot of code
        for ( my $i = 0 ; $i <= $#taxsz ; $i++ ) {
            last if $lastqs;
            $fam = $taxsz[$i];
            my $thissz = $sztaxa{$fam};
            if ( $thissz > $TNTB ) {
                if ( $thissz >= $TNTD ) {
                    print
    "family $fam has size $thissz - we think this is too big - skipping\n";
                    next;
                }
                if ( -s "data/$fam/oid.tre" ) {
                    print "family $fam alrady done\n";
                    next;
                }
                if ( $qsct >= $MAXQS ) {
                    $lastbig = $i;    #last entry used
                    print "max qsubs reached\n";
    
                    #                my @qskeys = keys %qsid;
                    open( SQ, ">$SCHEDLOG" );
                    print SQ "group $qsgrp\n";
    
                    #                foreach my $qs (@qskeys){
                    #                    print SQ "$qs\n" if ($qs != $qsgrp);
                    #                }
                    close SQ;
                    last;
                }
                if ( $actfam{$fam} ) {
                    print "skip $fam because active\n";
                    next;
                }
    
                #            if ($fam <= $bigest){
                #               my $fnd = 0;
                #               for my $famx(@actfam){
                #                   if ($fam == $famx){
                #                       print "skip $fam because active\n";
                #                       $fnd = 1;
                #                       last;
                #                   }
                #               }
                #               next if ($fnd); #cannot use this as already found
                #           }
    
                $qsct++ if ( startqs( \%qsid, 0, $fam ) );
                $lastbig++;
                $totqs++;
    
    #            my $walmem = setlimit($sztaxa{$fam});
    #            my $JOB_SCRIPT="$OID_USER_DIR/run_oid_job.sh";
    #            my $PARAMS = "-l nodes=1:ppn=$NCPU -l $walmem";
    #            my $PARAMS = "-l nodes=1:ppn=$NCPU -l mem=8GB,walltime=12:00:00";
    #            my $submit = q/qsub  /."$PARAMS $JOB_SCRIPT".q/ -v arg1="-at",arg2="^/."$fam".'$" 2>/dev/null';
    #            print "$submit\n";
    #    my $pid = open(QS,"qsub -l $PARAMS $JOB_SCRIPT -v arg1=/"$at/",arg2=/"^$fam.$dol/" 2>/dev/null |");
    #            my $pid = open(QS,"$submit |");
    #            my $qs = <QS>;
    #            if ($qs){
    #                chomp $qs;
    #                print "start qsid $qs for family $fam\n";
    #                $qsid{$qs} = $fam;
    #                $qsct++;
    #                $totqs++;
    #            } else {
    #                print "qsub failed\n";
    #            }
            }
            else {
                print "no more big jobs at $totqs\n";
                $lastqs = 1;
                last;
            }
        }
    ## now do serial stuff (later support multiple)
        my $famx;
        for $fam (@famtbl) {
            my $thissz = $sztaxa{$fam};
            if ( $thissz <= $TNTA ) {
                if ( $thissz < 4 ) {
                    print "skipping fam $fam too few taxa $thissz\n";
                    next;
                }
                if ( -s "data/$fam/oid.tre" ) {
                    print "$fam already done\n";
                    next;
                }
                my $famdol = '^' . "$fam" . '$';
                alignFamily("$famdol");
                makeTree("$famdol");
                chkuptm();
                ( $oldqs, $qsfam ) = testqs( \%qsid );
                if ($oldqs) {
                    print "completed family $qsfam qcnt $qsct\n";
                    $qsct--;
                    if ( $oldqs == $qsgrp and !-f "$PLDONE" ) {
                        my $qs = strtgrp( $ngroup, $grouplo, $grouphi );
                        if ($qs) {
                            $qsid{$qs} = 1;
                            $qsgrp = $qs;
                            $qsct++;
                            print "start qsid $qs for  group\n";
    
                            #                           my @qskeys = keys %qsid;
                            open( SQ, ">$SCHEDLOG" );
                            print SQ "group $qsgrp\n";
                            print SQ "bigmon $qsbigmon\n";
    
                       #                           foreach my $qs (@qskeys){
                       #                              print SQ "$qs\n" if ($qs != $qsgrp);
                       #                           }
                            close SQ;
                        }
                        else {
                            print "could not start group\n";
    
                            #                             $qsct--;
                        }
                        next;
                    }
                    if ( $oldqs == $qsbigmon and !-f "$BMDONE" ) {
                        my $qs = strbmon();
                        if ($qs) {
                            $qsid{$qs} = 1;
                            $qsbigmon = $qs;
                            $qsct++;
                            print "start qsid $qs for  Bigmon\n";
    
                            #                           my @qskeys = keys %qsid;
                            open( SQ, ">$SCHEDLOG" );
                            print SQ "group $qsgrp\n";
                            print SQ "bigmon $qsbigmon\n";
    
                       #                           foreach my $qs (@qskeys){
                       #                              print SQ "$qs\n" if ($qs != $qsgrp);
                       #                           }
                            close SQ;
                        }
                        else {
                            print "could not start Bigmon\n";
    
                            #                             $qsct--;
                        }
                        next;
                    }
                    while ( $qsct <= $MAXQS and !$lastqs ) {
    
                        $famx = $taxsz[$lastbig];
                        if ( $sztaxa{$famx} <= $TNTB ) {
                            print "no more big jobs at $lastbig\n";
                            $lastqs = 1;    ##no more
    
                            #                        $qsct--;
                            next;
                        }
                        if ( -s "data/$famx/oid.tre" ) {
                            print "family $famx alrady done\n";
                            $lastbig++;
                            next;
                        }
                        if ( $actfam{$famx} ) {
                            print "skip $famx because active\n";
                            next;
                        }
    
           #                    if ($famx<= $bigest){
           #                       my $fnd = 0;
           #                       for my $famy(@actfam){
           #                           if ($famx == $famy){
           #                               print "skip $famx because active\n";
           #                               $lastbig++;
           #                               $fnd = 1;
           #                               last;
           #                           }
           #                       }
           #                       next if ($fnd); #cannot use this as already found
           #                    }
                        $qsct++ if ( startqs( \%qsid, $oldqs, $famx ) );
                        $lastbig++;
    
                       #                    my @qskeys = keys %qsid;
                       #                    open (SQ, ">$SCHEDLOG");
                       #                    print SQ "group $qsgrp\n";
                       #                    foreach my $qs (@qskeys){
                       #                        print SQ "$qs\n" if ($qs != $qsgrp);
                       #                    }
                       #                    close SQ;
                        $totqs++;
                    }    ## for inner while
                }    ## if $oldqs
            }
            else {
                print "did all sequential align/maketre\n";
                last;
            }
        }    ## step forward through small taxa
        print "end of seq with qsct $qsct\n";
    
        #                my @qskeys = keys %qsid;
        open( SQ, ">$SCHEDLOG" );
        print SQ "group $qsgrp\n";
        print SQ "bigmon $qsbigmon\n";
    
        #                foreach my $qs (@qskeys){
        #                    print SQ "$qs\n" if ($qs != $qsgrp);
        #                }
        close SQ;
        sleep 15;           #give a little time for  things to settle
        %qsid   = ();       #empty
        $fndme  = 0;
        %actfam = ();       #empty
        until ($fndme) {    ## sometimes qstat doesn't work - so wont find me
            ( $fndme, %actfam ) = activqs( \%qsid, $qsgrp, $qsbigmon );
        }
    
        #        $bigest = 999999;
        #        if (@actfam){
        #            @actfam = sort {$a <=> $b} @actfam;
        #            $bigest = $actfam[$#actfam];
        #        }
        $qsct = keys %qsid;    #init with actual active tasks
        print "$qsct qsid \n";
        my $tmnw = localtime;
        print "$tmnw\n";
        $lastbig = 1;
    
        #    $lastqs = 0;  #spin through data 1 more time im case we missed one
        my $flag = 0;
        $lastqs = 1;           # must keep 1 because bigmon does all qsubs
        while ( $qsct > 0 ) {
            ( $oldqs, $qsfam ) = testqs( \%qsid );
            $flag = 0;
            if ($oldqs) {
                print "complete fam $qsfam qcnt $qsct\n";
                $qsct--;
                $flag = 1;
                if ( $oldqs == $qsgrp and !-f "$PLDONE" ) {
                    my $qs = strtgrp( $ngroup, $grouplo, $grouphi );
                    if ($qs) {
                        $qsid{$qs} = 1;
                        $qsct++;
                        $qsgrp = $qs;
                        print "start qsid $qs for  group\n";
    
                        #                    my @qskeys = keys %qsid;
                        open( SQ, ">$SCHEDLOG" );
                        print SQ "group $qsgrp\n";
    
                        #                    foreach my $qs (@qskeys){
                        #                       print SQ "$qs\n" if ($qs != $qsgrp);
                        #                    }
                        close SQ;
                    }
                    else {
                        print "could not start group\n";
    
                        #                     $qsct--;
                    }
                    next;
                }
                if ( $oldqs == $qsbigmon and !-f "$BMDONE" ) {
                    my $qs = strbmon();
                    if ($qs) {
                        $qsid{$qs} = 1;
                        $qsbigmon = $qs;
                        $qsct++;
                        print "start qsid $qs for  Bigmon\n";
                        open( SQ, ">$SCHEDLOG" );
                        print SQ "group $qsgrp\n";
                        print SQ "bigmon $qsbigmon\n";
                        close SQ;
                    }
                    else {
                        print "could not start Bigmon\n";
                    }
                    next;
                }
                while ( $qsct <= $MAXQS and !$lastqs ) {
                    $famx = $taxsz[$lastbig];
                    if ( $sztaxa{$famx} >= $TNTD ) {   # we skip the really big ones
                        $lastbig++;
                        next;
                    }
                    if ( $sztaxa{$famx} <= $TNTB ) {
                        print "no more big jobs at $lastbig\n";
                        $lastqs = 1;                   ##no more
                        last;
                    }
                    if ( -s "data/$famx/oid.tre" ) {
                        print "family $famx alrady done\n";
                        $lastbig++;
                        next;
                    }
                    if ( $actfam{$famx} ) {
                        print "skip $famx because active\n";
                        next;
                    }
    
           #                    if ($famx<= $bigest){
           #                       my $fnd = 0;
           #                       for my $famy(@actfam){
           #                           if ($famx == $famy){
           #                               print "skip $famx because active\n";
           #                               $lastbig++;
           #                               $fnd = 1;
           #                               last;
           #                           }
           #                       }
           #                       next if ($fnd); #cannot use this as already found
           #                    }
                    $qsct++ if ( startqs( \%qsid, $oldqs, $famx ) );
                    $lastbig++;
                    $totqs++;
                    next;
                }    # while
            }
            next if ($flag);
    
            ## later put in code to restart outselves if we are too old
            print "sleeping\n";
            sleep 300;
            chkuptm();
    
            #        $tmnw = time();
            #        my $td = $tmnw - $tstrt;
            #        if ($td > 72000) {
            #            print "time now ",localtime($tmnw) ," diff ",$td ,"\n";
            #            open (SQ,">schedrst");
            #            print SQ "$lastbig";
            #            close(SQ);
            #            restartme();  #restarts my task
            #            }
    
        }
    
        open( SQ, ">$SCHEDLOG" );
        print SQ "group $qsgrp\n";
        print SQ "bigmon $qsbigmon\n";
        close SQ;
        print "all $totqs qsubs ended\n";
        open( FI, ">$SCHEDONE" );
        print FI "$totqs\n";
        close FI;
    }
    exit(0);
}    #end sub

#
# Find all ortholog sets
# Given a parenthetical tree, return list of ortholog groups (list of list refs).
#
sub findOrthologs($) {
    my $dirRE = shift;
    $dirRE = '.*' if !defined($dirRE);
    my $verbose      = $VERBOSE;
    my $thisFunc     = ( caller(0) )[3];
    my $treeFile     = "oid.tre";
    my $orthologFile = "orthologs";

    # Change to data dir
    my $savDir = getcwd;
    chdir $OID_DATADIR;

    # Go over each family
    my $count = 0;
    print "Generating ortholog groups ...\n" if $verbose;
    foreach my $dir (<[1-9]*>) {
        next if $dir !~ /$dirRE/;
        $count++;
        if ( $verbose == 1 ) {
            $| = 1;    # flush STDOUT
            print "\b" x 7;
            printf "%6d", $count;
            $| = 0;
        }
        elsif ( $verbose == 2 ) {
            print "Finding orthologs for family $dir\n";
        }

        chdir $dir
          or die "$thisFunc: failed to change to family directory $dir: $!\n";

        if ( !-f $treeFile || -f $orthologFile ) {
            chdir $OID_DATADIR;
            next;
        }

        #my $parenTree = (untranslateTree($treeFile))[0];
        my $parenTree = ( parentheticalTree($treeFile) )[0];
        my $phyTree =
          new OrthologID::PhyloTree( $parenTree, [ ( @INGROUP, @OUTGROUP ) ] );
        my @orthGroups = $phyTree->orthologGroups( @INGROUP, @OUTGROUP );
        open FH, ">$orthologFile"
          or die "$thisFunc: failed to open ortholog file for writing: $!\n";
        foreach my $gp (@orthGroups) {
            my @g = @$gp;
            print FH "@g\n";
        }
        close FH;
        chdir $OID_DATADIR;

    }
    print "\n" if $verbose == 1;
    chdir $savDir;
}
1;

__END__


=head1 NAME

OrthologID - Collection of subroutines for the OrthologID framework

=head1 SYNOPSIS

  use OrthologID;


=head1 DESCRIPTION

OrthologID is a parsimony-based ortholog identification framework.
The OrthologID module contains all essential subroutines for
ortholog identification, including all-against-all blast, clustering,
alignment, tree searching, and ortholog extractions.  These subroutines are
intended to be used as part of a pipeline, as in the OrthologID software
package.


=head1 SEE ALSO

Chiu, J. C., Lee, E. K., Egan, M. G., Sarkar, I. N., Coruzzi, G. M., 
and Desalle, R. 2006. OrthologID: automation of genome-scale ortholog 
identification within a parsimony framework. Bioinformatics 22, 6 
(Mar. 2006), 699-707. DOI= http://dx.doi.org/10.1093/bioinformatics/btk040 

=head1 AUTHOR

Ernest K Lee, E<lt>elee@amnh.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2007-2010 by American Museum of Natural History

=cut
