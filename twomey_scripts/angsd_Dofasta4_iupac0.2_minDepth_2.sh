for file in *.bam
do
    angsd -doFasta 4 -doCounts 1 -minQ 20 -minMapQ 30 -remove_bads 1 -uniqueOnly -setMinDepth 2 -iupacRatio 0.2 -i $file
    gunzip angsdput.fa.gz
    sed -i "s/uce/$file|uce/g" angsdput.fa
    mv angsdput.fa $file-uce-angsd-dofasta.fasta
    rm angsdput.arg
done
