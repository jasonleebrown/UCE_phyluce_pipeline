for dir in */split-adapter-quality-trimmed; do
    cd $dir
    cd ..
    sample=`basename "$(pwd)"`
    print=$sample
    cd split-adapter-quality-trimmed
    bwa mem /home/bender/Desktop/wilson_bassleri_project_042820/mtGenome/reference/mtgenome_reference3.fasta *READ1.fastq.gz *READ2.fastq.gz -B 2 -t 6 > $sample.sam
    samtools view -S -u $sample.sam -o $sample.bam
    ##if this errors out use following code
    samtools sort $sample.bam -o $sample.sorted.bam
    ##samtools sort $sample.bam $sample.sorted
    rm $sample.sam
    rm $sample.bam
    mv $sample.sorted.bam /home/bender/Desktop/wilson_bassleri_project_042820/mtGenome/angsd_bams/$sample.bam
    cd ../../
done