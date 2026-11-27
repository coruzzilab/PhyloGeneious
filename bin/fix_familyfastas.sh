#!/bin/sh

grep -q ".>" $1
#grep ">" | grep -q -v "^>"

if [[ $?==0 ]]; then
sed -i "s/\(.\)>/\1\n>/g;s/\f//" $1;
fi
