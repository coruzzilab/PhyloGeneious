shopt -s extglob 

TEST_DIR=/scratch/$USER/testdata/Q20_all
cd $TEST_DIR/data

find . ! -name 'FAMILY' -type f -exec rm -f {} +
rm -f $TEST_DIR/log/job/bigmdone
touch .family.done