#!/usr/bin/env perl
#
#
# Copyright (C) 2006-2011 Ernest K. Lee
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
# Call different OrthologID function according to given options and argument.
#

use warnings;

my $OID_HOME;
BEGIN {
	$OID_HOME=$ENV{'OID_HOME'};
	die "Environment variable OID_HOME is not defined ... exiting.\n" if ! defined($OID_HOME);
	$OID_USER_DIR = $ENV{'OID_USER_DIR'};
	die "Environment variable OID_USER_DIR is not defined ... exiting.\n"
	  if !defined($OID_USER_DIR);
}

use lib "$OID_HOME/lib";
use OrthologID; 
use Getopt::Std;
use strict;
print "this program checks for tnt errors and tries to restart tnt\n";
my %locproc  = ();
my @lcutoff = ();
my $prkey;
#package OrthlogID;
makeTree(3);
print "TNTA $TNTA\n";
foreach $prkey(keys %{procptr}){
    print "proc key $prkey\n";
    $locproc{$prkey} = ${procptr}{$prkey};
}
#@lcutoff = @{OrthologID::cutoff};
my $nkeys = keys %locproc;
foreach my $val (@{OrthologID::cutoff}){
    push @lcutoff,$val;
    print "cutoff val $val\n";
}
print "there are $nkeys in local procptr\n";
#print "cutoff" ,join ",",@lcutoff,"\n";

