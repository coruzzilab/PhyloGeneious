#!/bin/bash
cd $OID_USER_DIR
$OID_HOME/bin/topshell.pl
chmod +x pipe.sh
$OID_USER_DIR/pipe.sh $@
