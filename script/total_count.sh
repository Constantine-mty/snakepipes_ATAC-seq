
for bam in ../bam/*sortp_genome_rmchrM_rmdup.bam; do
    sample=$(basename $bam .sortp_genome_rmchrM_rmdup.bam)
    samtools idxstats $bam > ../count_table/${sample}.total.count_table
done