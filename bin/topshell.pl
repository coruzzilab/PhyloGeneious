#!/usr/bin/env perl
use strict;
my $OID_HOME;
my $OID_USER_DIR;
my $OID_WRAPPER;

BEGIN {
    $OID_HOME = $ENV{'OID_HOME'};
    die "Environment variable OID_HOME is not defined ... exiting.\n"
      if !defined($OID_HOME);
    $OID_USER_DIR = $ENV{'OID_USER_DIR'};
    die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
      if !defined($OID_USER_DIR);
    $OID_WRAPPER = $ENV{'OID_WRAPPER'};
}

use lib "$OID_HOME/lib";
use OrthologID;
my $HPC;
$HPC = $OrthologID::HPC;      #env nogood
$HPC = 'P' if ( not $HPC );
print "we are in hpc $HPC\n";
open FO, ">pipe.sh" or die "cannot create pipe.sh";
print FO "#!/bin/bash\n";
print FO "export HPC=$HPC\n";
print FO "if [[ ! -d toplog ]]; then mkdir toplog; fi\n";

if ( $HPC =~ /S/ ) {
    print FO "sbatch -o toplog/%J.out \$@ $OID_WRAPPER $OID_HOME/bin/pipe.slu\n";
}
elsif ( $HPC =~ /P/ ) {
    print FO "qsub \$@ $OID_WRAPPER $OID_HOME/bin/pipe.pbs\n";
}
else {
    print FO "qsub \$@ $OID_WRAPPER $OID_HOME/bin/pipe.pbs\n";
}
print "created pipe.sh\n";
close FO;
