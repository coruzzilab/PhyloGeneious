#use: make_tfsubproc.sh /path/to/FAMILY.aligned.revfasta k RAND

TESTDATA=/scratch/$USER/testdata

BUFFER=$'log +pi'
BUFFER+=$2
BUFFER+=$'.log\nmxram 15000;\ntaxname=;\ntaxname +64;\nnstates DNA;\nwatch=;\np &'
BUFFER+=$1
BUFFER+=$';\nhold 11;\nrs '
BUFFER+=$3
BUFFER+=$';\ncol 3;\npfijo: slack 100 numsels 1;\nxi;\nriddup;\nmu=rep1ho10wag;\nriddup-;\ntsave pi'
BUFFER+=$2
BUFFER+=$'.rat;\nsave.;\ntimeout 4:00:00;\nhold 10;\nquote ------------------ pinion (pfijo) SEARCH 1 of 1;\nxi;\nriddup;\npfijo 0;\nriddup-;\nsave.;\ntsave/;\nk0;\nhold 40;\nshort pi' 
BUFFER+=$2
BUFFER+=$'.rat;\nbest;\nexport - pi'
BUFFER+=$2
BUFFER+=$'.tre;\nlog/;\nzz;\nproc/;'

echo "$BUFFER" > oblong_pfijo${2}.proc
