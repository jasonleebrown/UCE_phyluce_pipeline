for file in *.bam
do
    angsd -doFasta 2 -doCounts 1 -minQ 10 -minMapQ 20 -remove_bads 1 -uniqueOnly -setMinDepth 0.5 -i $file
    gunzip angsdput.fa.gz
	##change "T627178.1" to match the file name at top of output fasta files
    sed -i "s/MT627178_AF2673_Ameerega_hahneli/$file|mtDNA/g" angsdput.fa
    mv angsdput.fa $file-angsd.fasta
    rm angsdput.arg
done
