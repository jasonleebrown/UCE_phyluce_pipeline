#!/bin/bash

#finds the number of UCEs in the file by grepping and wc'ing the number of "charset" lines
uce_num=`grep charset $1 | wc -l`
echo "There are "$uce_num" UCEs"

#make an output file that is a copy of the input file
cp $1 correctedpartitions.tmp

#change - dashes in UCE names to _ underscores
sed -ri "s/(charset 'uce)-(.*;)/\1_\2/g" correctedpartitions.tmp

#make a file consisting of only "charset" lines
grep charset correctedpartitions.tmp > tmp1.tmp

#edit the file with sed so that it only consists of the relevant UCE loci
sed -ri "s/charset '(uce_.*)' = .*;/\1/g" tmp1.tmp

#start the charpartition construct
echo "charpartition loci = " > tmp2.tmp

#initialize an index variable, which we will increment in the following for loop
index=1

#for loop writes to tmp2tmp
for i in $(cat tmp1.tmp); do echo "$index:$i, " >> tmp2.tmp; index=`expr $index + 1`; done

#use tr to delete all newlines in tmp2.tmp and make tmp3.txt
cat tmp2.tmp | tr -d "\n" > tmp3.tmp

#replace final comma in charpartition.txt with a semicolon and END; construct on the next line
sed -ri 's/(.*),/\1;\nEND;/g' tmp3.tmp

#use head to get rid of previous (wrong) charpartition construct and create new file
head -n -2 correctedpartitions.tmp > correctedpartitions.nexus

#append contents of charpartition.txt to the output .nexus file
cat tmp3.tmp >> correctedpartitions.nexus

#remove single quotes from the output file so that charsets match charpartitions
sed -i "s/'//g" correctedpartitions.nexus

echo "job complete; see correctedpartitions.nexus"

rm tmp1.tmp
rm tmp2.tmp
rm tmp3.tmp
rm correctedpartitions.tmp
