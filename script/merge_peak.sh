
cat ../result/peak_calling/*overlapped_peaks.narrowPeak | grep -v chrY | grep -v chrM > ../result/downstream/diff_peak/merge_atac.peak

# 排序
bedtools sort -i ../result/downstream/diff_peak/merge_atac.peak -g ../ref_data/mm10/mm10_no_alt_analysis_set_ENCODE.fasta.fai > ../result/downstream/diff_peak/merge_atac.sort.peak

# 查看排序
cut -f 1 ../result/downstream/diff_peak/merge_atac.sort.peak | uniq

# 默认的bedtools merge
bedtools merge -i ../result/downstream/diff_peak/merge_atac.sort.peak > ../result/downstream/diff_peak/merge_atac.sort.combine.peak

# 对合并后的Peak进行定量，每个bam重复计算
for bam in ../result/bam/*sortp_genome_rmchrM_rmdup.bam; do
    # 取文件名（不带路径和后缀）
    fname=$(basename "$bam")
    sample=${fname%.sortp_genome_rmchrM_rmdup.bam}

    bedtools coverage \
        -a ../result/downstream/diff_peak/merge_atac.sort.combine.peak \
        -b "$bam" \
        -counts \
        -g ../ref_data/mm10/mm10_no_alt_analysis_set_ENCODE.fasta.fai \
        -sorted \
        > "../result/downstream/diff_peak/${sample}_combine_peak.count.table"
done

# 检查count table文件行数是否一致
wc -l ../result/downstream/diff_peak/*_combine_peak.count.table


