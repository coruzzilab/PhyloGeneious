#!/bin/bash
if [ ! -f $OID_HOME/.tntagreed ]; then
echo "ERROR: PhyloGeneious not set up"
echo "Please run the setup_license.sh script in the pipeline installation directory before proceeding"
exit 1
fi

cd $OID_USER_DIR
$OID_HOME/bin/topshell.pl #$ENV_WRAPPER 
chmod +x pipe.sh
$OID_USER_DIR/pipe.sh $@
