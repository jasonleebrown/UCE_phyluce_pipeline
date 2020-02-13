#!/bin/bash

#grab all character blocks and put them in a temporary file
grep uce $1 > tmp1.tmp

#add charset qualifier to each line
sed -i 's/uce/charset uce/g' tmp1.tmp

#initialize charpartition temp file
echo "charpartition swscen = " > tmp2.tmp

#rearrange charsets and append to tmp2.tmp
#sed -r 's/charset (uce.*) = (.*);/\1: \2, /g' tmp1.tmp >> tmp2.tmp
sed -r 's/charset (uce.*) = (.*);/GTR:\1, /g' tmp1.tmp >> tmp2.tmp

#use tr to remove newlines from tmp2.tmp and move to tmp3.tmp
cat tmp2.tmp | tr -d "\n" > tmp3.tmp

#use sed to replace final comma with semicolon and END; construct
sed -ri 's/(.*),/\1;\nEND;/g' tmp3.tmp

#make final output file
echo "#nexus" > swscen_for_iqtree.nexus
echo "begin sets;" >> swscen_for_iqtree.nexus
cat tmp1.tmp >> swscen_for_iqtree.nexus
cat tmp3.tmp >> swscen_for_iqtree.nexus

rm tmp1.tmp
rm tmp2.tmp
rm tmp3.tmp
