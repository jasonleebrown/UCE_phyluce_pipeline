mkdir consensus
for file in *.fasta
do
    filename=`basename $file .fasta`
    cons -sequence $file -outseq consensus/$filename-con.fasta -name $filename-con
done