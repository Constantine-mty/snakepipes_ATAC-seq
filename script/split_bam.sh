
#mkdir ../result/split_bam

#for bam in ../bam/*sortp_genome_rmchrM_rmdup.bam; do
#    # select different length
#    samtools view -@ 5 -h $bam | \
#       awk 'substr($0,1,1)=="@" || ($9>= 120) || ($9<=-120)' | \
#       samtools view -@ 5 -hb > split_bam/K562-ATACSeq-rep1.ENCFF534DCE.Over120bp.bam &
#
#
#    samtools view -@ 5 -h $bam | \
#       awk 'substr($0,1,1)=="@" || ($9 < 120 && $9 > 0) || ($9 > -120 && $9 < 0)' | \
#       samtools view -@ 5 -hb > split_bam/K562-ATACSeq-rep1.ENCFF534DCE.Shorter120bp.bam &
#done


# select different length
samtools view -@ 5 -h ../result/bam/Act_FF_1.sortp_genome_rmchrM_rmdup.bam | \
    awk 'substr($0,1,1)=="@" || ($9>= 120) || ($9<=-120)' | \
    samtools view -@ 5 -hb > ../result/split_bam/Act_FF_1.Over120bp.bam &


samtools view -@ 5 -h ../result/bam/Act_FF_1.sortp_genome_rmchrM_rmdup.bam | \
    awk 'substr($0,1,1)=="@" || ($9 < 120 && $9 > 0) || ($9 > -120 && $9 < 0)' | \
    samtools view -@ 5 -hb > ../result/split_bam/Act_FF_1.Shorter120bp.bam &



