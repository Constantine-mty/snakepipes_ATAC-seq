
rgt-hint footprinting --atac-seq --paired-end --output-prefix=FF1 --organism=mm10 \
    /media/desk16/MaTianyu/Pipeline/snakepipes_ATAC-seq/result/bam/Act_FF_1.sortp_genome_rmchrM_rmdup.bam \
    /media/desk16/MaTianyu/Pipeline/snakepipes_ATAC-seq/result/peak_calling/Act_FF_1.genrich_peaks.rmBlacklist.narrowPeak


for celltype in FF1 FF2 FF3 KO1 KO2 KO3;
do
    rgt-hint footprinting --atac-seq --paired-end --organism=mm10 \
    --output-location=./downstream/footprints \
    --output-prefix=${celltype} \
    ./BAM/${celltype}.bam ./Peaks/${celltype}_peaks.narrowPeak
done


# Footprints region Motif scan

# mkdir MotifMatching

# for celltype in CLP CMP GMP HSC LMPP MEP MPP pDC;
# do
# rgt-motifanalysis matching \
# --organism=hg19 \
# --output-location=./MotifMatching \
# --input-files ./Footprints/${celltype}.bed
# done


# Differential analysis footprint signal

# mkdir -p DiffFootprinting

# rgt-hint differential --organism=hg19 --bc --nc 64 \
# --mpbs-files=./MotifMatching/CLP_mpbs.bed,\
# ./MotifMatching/CMP_mpbs.bed,\
# ./MotifMatching/GMP_mpbs.bed,\
# ./MotifMatching/HSC_mpbs.bed,\
# ./MotifMatching/LMPP_mpbs.bed,\
# ./MotifMatching/MEP_mpbs.bed,\
# ./MotifMatching/MPP_mpbs.bed,\
# ./MotifMatching/pDC_mpbs.bed \
# --reads-files=./BAM/CLP.bam,./BAM/CMP.bam,./BAM/GMP.bam,\
# ./BAM/HSC.bam,./BAM/LMPP.bam,./BAM/MEP.bam,./BAM/MPP.bam,\
# ./BAM/pDC.bam \
# --conditions=CLP,CMP,GMP,HSC,LMPP,MEP,MPP,pDC \
# --output-location=./DiffFootprinting --output-prefix=All