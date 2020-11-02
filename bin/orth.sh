date
time
# Extract orthologs
echo "Extracting orthologs ..."
$OID_BIN/orthologid.pl -O
date
time

# Generate big matrix
echo 'Generating matrix ...'
$OID_BIN/orth2matrix.pl

date
time
