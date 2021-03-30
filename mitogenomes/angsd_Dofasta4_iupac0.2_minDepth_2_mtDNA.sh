for file in *.bam
do
    angsd -doFasta 2 -doCounts 1 -minQ 20 -minMapQ 30 -remove_bads 1 -uniqueOnly -setMinDepth 1 -i $file
    gunzip angsdput.fa.gz
	##change "T627178.1" to match the file name at top of output fasta files
    sed -i "s/T627178.1/$file|mtDNA/g" angsdput.fa
    mv angsdput.fa $file-angsd.fasta
    rm angsdput.arg
done