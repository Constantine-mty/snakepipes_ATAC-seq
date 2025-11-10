#!/bin/bash

bam_dir="../bam"   # 根据实际路径修改
bam_files=($bam_dir/*sortp_genome_rmchrM_rmdup.bam)

# 检查是否至少有两个文件
if [ ${#bam_files[@]} -lt 2 ]; then
  echo "Need at least 2 BAM files for multiBamSummary!"
  exit 1
fi

# 生成 labels
labels=()
for f in "${bam_files[@]}"; do
    base=$(basename "$f")
    label=${base%%.sortp_genome_rmchrM_rmdup.bam}
    labels+=("$label")
done

# 拼接成字符串
bamfiles_str=$(IFS=$' '; echo "${bam_files[*]}")
labels_str=$(IFS=$' '; echo "${labels[*]}")

# 输出命令（可直接运行）

plotFingerprint \
 --bamfiles $bamfiles_str \
 --binSize=1000 \
 --plotFile ATAC.fingerprint.pdf \
 --labels $labels_str \
 -p 8 \
 &> fingerprint.log

# multiBamSummary bins \
#  --bamfiles $bamfiles_str \
#  --labels $labels_str \
#  --outFileName multiBamArray_ATAC.npz \
#  --binSize 5000 \
#  -p 8 \
#  &> multiBamSummary.log

# plotCorrelation \
#  --corData multiBamArray_ATAC.npz \
#  --plotFile ATAC_correlation_bin.pdf \
#  --outFileCorMatrix ATAC_correlation_bin.txt \
#  --whatToPlot heatmap \
#  --corMethod spearman

