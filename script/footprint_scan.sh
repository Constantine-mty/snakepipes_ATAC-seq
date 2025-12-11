
#rgt-hint footprinting --atac-seq --paired-end --output-prefix=FF1 --organism=mm10 \
#    /media/desk16/MaTianyu/Pipeline/snakepipes_ATAC-seq/result/bam/Act_FF_1.sortp_genome_rmchrM_rmdup.bam \
#    /media/desk16/MaTianyu/Pipeline/snakepipes_ATAC-seq/result/peak_calling/Act_FF_1.genrich_peaks.rmBlacklist.narrowPeak


# for celltype in Act_FF_1 Act_FF_2 Act_FF_3 Act_KO_1 Act_KO_2 Act_KO_3;
# do
#     rgt-hint footprinting --atac-seq --paired-end --organism=mm10 \
#     --output-location=./downstream/footprint \
#     --output-prefix=${celltype} \
#     ./bam/${celltype}.sortp_genome_rmchrM_rmdup.bam ./peak_calling/${celltype}.macs2_peaks.rmBlacklist.narrowPeak
# done


# Footprints region Motif scan

# mkdir MotifMatching

# for celltype in Act_FF_1 Act_FF_2 Act_FF_3 Act_KO_1 Act_KO_2 Act_KO_3;
# do
#     rgt-motifanalysis matching \
#     --organism=mm10 \
#     --output-location=./MotifMatching \
#     --input-files ./${celltype}.bed
# done


# Differential analysis footprint signal

#mkdir -p DiffFootprinting

rgt-hint differential --organism=hg19 --bc --nc 64 \
--mpbs-files=./MotifMatching/Act_FF_1_mpbs.bed,\
./MotifMatching/Act_FF_2_mpbs.bed,\
./MotifMatching/Act_FF_3_mpbs.bed,\
./MotifMatching/Act_KO_1_mpbs.bed,\
./MotifMatching/Act_KO_2_mpbs.bed,\
./MotifMatching/Act_KO_3_mpbs.bed \
--reads-files=../../bam/Act_FF_1.sortp_genome_rmchrM_rmdup.bam,\
../../bam/Act_FF_2.sortp_genome_rmchrM_rmdup.bam,\
../../bam/Act_FF_3.sortp_genome_rmchrM_rmdup.bam,\
../../bam/Act_KO_1.sortp_genome_rmchrM_rmdup.bam,\
../../bam/Act_KO_2.sortp_genome_rmchrM_rmdup.bam,\
../../bam/Act_KO_3.sortp_genome_rmchrM_rmdup.bam \
--conditions=Act_FF_1,Act_FF_2,Act_FF_3,Act_KO_1,Act_KO_2,Act_KO_3 \
--output-location=./DiffFootprinting --output-prefix=All