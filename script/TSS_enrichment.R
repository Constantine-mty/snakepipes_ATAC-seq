library(readr)
library(dplyr)
library(stringr)

# -------------------------------
# 参数：count_table目录和物种
# -------------------------------
count_dir <- "../count_table/"
peak_dir <- "../peak_calling"
species <- "mouse"  # "human" 或 "mouse"

# 根据物种设置参考基因组大小和染色体列表
if(species == "mouse"){
  genome_size <- 2.7e9
  chr_list <- paste0("chr", c(1:19, "X", "Y"))
} else if(species == "human"){
  genome_size <- 3.1e9
  chr_list <- paste0("chr", c(1:22, "X", "Y"))
} else{
  stop("Unknown species!")
}

# -------------------------------
# 自动获取样本列表
# -------------------------------
tss_files <- list.files(count_dir, pattern = "\\.TSS_2Kbp\\.count_table$", full.names = TRUE)
samples <- str_replace(basename(tss_files), "\\.TSS_2Kbp\\.count_table$", "")

# -------------------------------
# 初始化结果表
# -------------------------------

# 初始化结果表
results <- data.frame(
  Sample = samples,
  FRiP_peaks = NA,
  FRiP_summary = NA,
  TSS_enrichment = NA
)

# -------------------------------
# 循环处理每个样本
# -------------------------------


for (i in seq_along(samples)) {
  s <- samples[i]
  
  ## --------------------------
  ## FRiP 方法1: 原始 peak counts
  ## --------------------------
  peak_file <- file.path(peak_dir, paste0(s, ".overlapped_peaks.counts"))
  total_count_file <- file.path(count_dir, paste0(s, ".total.count_table"))
  
  if (file.exists(peak_file) & file.exists(total_count_file)) {
    peak_count <- read_tsv(peak_file, col_names = TRUE, skip = 1)
    colnames(peak_count) <- c("PeakID","Chr", "Start", "End", "Strand", "Length", "Count")
    total_peak_count <- sum(peak_count$Count)
    
    total_count_df <- read_tsv(total_count_file, col_names = F)
    colnames(total_count_df) <- c("chrom", "chr_len", "mapped_count", "unmapped_count")
    total_count <- sum(total_count_df$mapped_count)
    
    results$FRiP_peaks[i] <- total_peak_count / total_count
  }
  
  ## --------------------------
  ## FRiP 方法2: summary 文件
  ## --------------------------
  summary_file <- file.path(peak_dir, paste0(s, ".overlapped_peaks.counts.summary"))
  if (file.exists(summary_file)) {
    count_summary <- read_tsv(summary_file, col_names = TRUE)
    colnames(count_summary) <- c("AlignCategory", "Number")
    # 注意 Assigned 拼写
    results$FRiP_summary[i] <- count_summary$Number[count_summary$AlignCategory == "Assigned"] / sum(count_summary$Number)
  }
  
  ## --------------------------
  ## TSS enrichment
  ## --------------------------
  tss_file <- file.path(count_dir, paste0(s, ".TSS_2Kbp.count_table"))
  tss_2kb_count_df <- read_tsv(tss_file, col_names = F)
  colnames(tss_2kb_count_df) <- c("chrom", "tss_up_2k", "tss_down_2k", "gene_id", "tss", "strand", "region_count")
  
  # 根据物种选择 chr_list 和 genome_size
  genome_size <- 2.7e9  # 人
  chr_list <- paste0("chr", c(1:22, "X", "Y"))
  # 示例：小鼠样本名开头为 mm 或者手动指定 species
  if (grepl("^mm", s)) {
    genome_size <- 2.7e9  # 可替换为 mm10 实际基因组大小
    chr_list <- paste0("chr", c(1:19, "X", "Y"))
  }
  
  total_count_df.filter <- filter(total_count_df, chrom %in% chr_list)
  genome_coverage <- sum(total_count_df.filter$mapped_count) * 150 / genome_size
  tss_coverage <- sum(tss_2kb_count_df$region_count) * 150 / (4000 * nrow(tss_2kb_count_df))
  results$TSS_enrichment[i] <- tss_coverage / genome_coverage
}

# -------------------------------
# 保存结果到文件
# -------------------------------

str(results)

write_tsv(results, "ATAC_QC_summary.tsv")
print(results)



# test for single sample -------------------------------------------------------


# -------------------------------------------------------------------->>>>>>>>
# Fraction of reads in peaks (FriP)
# -------------------------------------------------------------------->>>>>>>>

# peak_count <- read_tsv("../peak_calling/Act_KO_1.overlapped_peaks.counts", col_names = T, skip = 1)
# colnames(peak_count) <- c("PeakID","Chr", "Start", "End", "Strand", "Length", "Count")
# peak_count
# 
# 
# total_peak_count <- sum(peak_count$Count)
# 
# 
# total_count_df <- read_tsv(file = "../count_table/Act_KO_1.total.count_table", col_names = F)
# colnames(total_count_df) <- c("chrom", "chr_len", "mapped_count", "unmapped_count")
# 
# total_count <- sum(total_count_df$mapped_count)
# 
# FriP <- total_peak_count/total_count


# -------------------------------------------------------------------->>>>>>>>
# Fraction of reads in peaks (FriP) Method2
# -------------------------------------------------------------------->>>>>>>>

# count_summary <- read_tsv("../peak_calling/Act_KO_1.overlapped_peaks.counts.summary", col_names = T)
# colnames(count_summary) <- c("AlignCategory", "Number")
# count_summary
# 
# FRiP <- count_summary$Number[count_summary$AlignCategory == "Assigned"] / sum(count_summary$Number)
# 
# FRiP

# -------------------------------------------------------------------->>>>>>>>
# load TSS count and calculate ATAC-seq TSS enrichment score
# -------------------------------------------------------------------->>>>>>>>
# load table 
# total_count_df <- read_tsv(file = "../count_table/Act_KO_1.total.count_table", col_names = F)
# colnames(total_count_df) <- c("chrom", "chr_len", "mapped_count", "unmapped_count")
# 
# total_count_df.filter <- filter(
#   total_count_df,
#   chrom %in% chr_list
# )
# 
# genome_coverage = sum(total_count_df.filter$mapped_count) * 150 / 2.7e9
# genome_coverage
# 
# # TSS count
# tss_2kb_count_df <- read_tsv(file = "../count_table/Act_KO_1.TSS_2Kbp.count_table", col_names = F)
# colnames(tss_2kb_count_df) <- c("chrom", "tss_up_2k", "tss_down_2k", "gene_id", "tss", "strand", "region_count")
# 
# tss_coverage = sum(tss_2kb_count_df$region_count) * 150 / (4000 * nrow(tss_2kb_count_df)) # 150 是 测序的 read 长度（bp）
# 
# # TSS enrichment score
# tss_enrichment_score = tss_coverage / genome_coverage
# tss_enrichment_score




