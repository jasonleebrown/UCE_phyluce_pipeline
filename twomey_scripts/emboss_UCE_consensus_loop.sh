mkdir consensus
for file in *.nexus
do
    filename=`basename $file .nexus`
    cons -sequence $file -outseq consensus/$filename-con.fasta -name $filename-con
done
