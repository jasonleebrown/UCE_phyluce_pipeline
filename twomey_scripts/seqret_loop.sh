mkdir seqret_out
for file in *.fa
do
    seqret -auto -sequence $file -snucleotide1 -sformat1 pearson -osformat2 msf -stdout > $file.msf
    seqret -auto -sequence $file.msf -snucleotide1 -sformat1 msf -osformat2 fasta -stdout > seqret_out/$file.fasta
    rm *.msf
done


