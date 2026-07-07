#!/bin/sh

cat $1 | grep -q ".>"
#grep ">" | grep -q -v "^>"

if [[ $?==0 ]]; then
sed -i "s/\(.\)>/\1\n>/g" $1;
fi
