#use: make_tfproc.sh /path/to/FAMILY.aligned.revfasta /path/to/oi.rat /path/to/outdir
#ratbld6 use: make_tfproc.sh /path/to/FAMILY.aligned.revfasta 15G /path/to/outdir #_ratfiles

if [[ $2 == '15G' ]]; then
  BUFFER=$'log +pi1.log\nmxram 15000;'
  else BUFFER=$'log +oi1.log\nmxram 5000;'
fi
BUFFER+=$'\ntaxname=;\ntaxname +64;\nnstates DNA;\nwatch=;\nreport +180/10/10;\np &'
BUFFER+=$1
BUFFER+=$';\nhold 2;\nrs 1;'
#;\nset chk isfile '
#BUFFER+=$2
#BUFFER+=$';\np ' 
#if (\'chk\'==1) 
#BUFFER+=$2

OUT_DIR=$3 #;end
if [[ $2 == '15G' ]]; then
  BUFFER+=$'\nmacro [200;\nmacro=;\nvar: chk;\n'
  for k in $(seq 0 1 $4); do
        BUFFER+=$'set chk isfile '
        BUFFER+=$OUT_DIR
        BUFFER+=$'/pi'
        BUFFER+=$k
        BUFFER+=$'.rat;\n'
        BUFFER+=$'if ("chk"==1) p '
        BUFFER+=$OUT_DIR
        BUFFER+=$'/pi'
        BUFFER+=$k
        BUFFER+=$'.rat;end;\n'
  done
  BUFFER+=$'set chk ntrees+10;\nhold  \'chk\';\nmacro-;\ntf: rounds 1 min 7;'
  BUFFER+=$'\ntsave pi.rat;\ntimeout 08:00:00;\nsort;\nscore;\nxi;\nriddup;\ntf;\nriddup-;\nsave .;\nriddup;\ntf;\nriddup-;\nsave .;\nriddup;\ntf;\nriddup-;\nsave .;\nriddup;\ntf;\nriddup-;\nsave .;\nriddup;\ntf;\nriddup-;\nsave .;\ntsave/ ;\nsort;\nscore;\nbest;\nscore;\ntsave oids.tre;\nsave .;\ntsave/ ;\ntaxname ] ;\ntchoose 0 ;\nexport - '
else 
  BUFFER+=$'\np ' 
  BUFFER+=$2
  BUFFER+=$';\nmacro [200;\nmacro=;\nvar: chk;\nset chk ntrees+10;\nhold  \'chk\';\nmacro-;\ntf: rounds 1 min 5;'
  BUFFER+=$'\ntsave oi1.rat;\ntimeout 08:00:00;\nsort;\nscore;\nxi;\nriddup;\ntf;\nriddup-;\nsave .;\nriddup;\ntf;\nriddup-;\nsave .;\nriddup;\ntf;\nriddup-;\nsave .;\nriddup;\ntf;\nriddup-;\nsave .;\nriddup;\ntf;\nriddup-;\nsave .;\ntsave/ ;\nsort;\nscore;\nbest;\nscore;\ntsave oids.tre;\nsave .;\ntsave/ ;\ntaxname ] ;\ntchoose 0 ;\nexport - '
fi
BUFFER+=$'oid.tre;\nlog/ ;\nzzz;\nproc/;'

echo "$BUFFER" > $OUT_DIR/oblong_treefuse.proc
