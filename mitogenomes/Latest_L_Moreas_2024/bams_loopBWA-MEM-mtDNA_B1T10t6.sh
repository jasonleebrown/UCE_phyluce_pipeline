for dir in */split-adapter-quality-trimmed; do
    cd $dir
    cd ..
    sample=`basename "$(pwd)"`
    print=$sample
    cd split-adapter-quality-trimmed
    bwa mem /home/student/Desktop/TD/work_directory/reference_sets/amee.fasta *READ1.fastq.gz *READ2.fastq.gz -B 1 -T 10 -t 6 > $sample.sam
    samtools view -S -u $sample.sam -o $sample.bam
    ##if this errors out use following code
    samtools sort $sample.bam -o $sample.sorted.bam
    ##samtools sort $sample.bam $sample.sorted
    rm $sample.sam
    rm $sample.bam
    mv $sample.sorted.bam /home/student/Desktop/TD/work_directory/angsd_bams/$sample.bam
    cd ../../
done
