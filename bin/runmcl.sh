 sbatch -N 1 --ntasks-per-node 1 -c 1 -o toplog/%J.out -t 12:00:00 --mem 16GB $OID_HOME/bin/runmcl.pl 16GB 20
