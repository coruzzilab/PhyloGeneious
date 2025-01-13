#!/bin/sh

#nwk to tnt
FILE=$1
OUT=$2
for TREE in `sed "s/\[\&R\] //g" $FILE| grep "^("`; do
sed "s/\[\&R\] //g" $FILE| sed "s/[0-9]\.[0-9]*//g" | sed "s/[:]//g" | sed "s/,/ /g" > tmp; #get trees
i=0; 
for ID in `sed "s/\[\&R\] //g" $FILE | sed "s/[0-9]\.[0-9]*//g" | sed "s/[:();\*]//g" | sed "s/,/ /g"`; do #get sequence IDs
sed -i "s/${ID}/${i}/" tmp; 
i=$((i+1)); done;
cat tmp >> tmp2;
done

#print output
echo xread > $OUT
echo 0 $i >> $OUT
sed "s/\[\&R\] //g" $FILE| sed "s/[0-9]\.[0-9]*//g" | sed "s/[:()]//g" | sed "s/,/ /g" >> $OUT
echo "collapse none;" >> $OUT
echo "tread '$(grep -c "^(" tmp2) trees'">> $OUT
cat tmp2 >> $OUT
echo "proc/;" >> $OUT

#clean temp files
rm tmp2;
rm tmp;

#2023 Veronica Sondervan
