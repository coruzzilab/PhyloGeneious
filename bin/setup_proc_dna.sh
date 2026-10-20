#!/bin/bash

cd $OID_DATADIR/$1
perl $OID_HOME/bin/setup_proc_dna.pl $1
