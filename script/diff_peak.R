####################################################################################################
# Data 2025-11-10
# Author Ma Tianyu
# E-mail 647892189@qq.com
# ChIP-seq data analaysis
####################################################################################################
rm(list=ls())

setwd(dir="~/Pipeline/snakepipes_ATAC-seq/")

library(tidyverse)

# 读取Peak count table ----------------------------------
count_df.ff_rep1 <- read_tsv(file = "./result/downstream/diff_peak/Act_FF_1_combine_peak.count.table", col_names = F)
count_df.ff_rep2 <- read_tsv(file = "./result/downstream/diff_peak/Act_FF_2_combine_peak.count.table", col_names = F)
count_df.ff_rep3 <- read_tsv(file = "./result/downstream/diff_peak/Act_FF_3_combine_peak.count.table", col_names = F)

count_df.ko_rep1 <- read_tsv(file = "./result/downstream/diff_peak/Act_KO_1_combine_peak.count.table", col_names = F)
count_df.ko_rep2 <- read_tsv(file = "./result/downstream/diff_peak/Act_KO_2_combine_peak.count.table", col_names = F)
count_df.ko_rep3 <- read_tsv(file = "./result/downstream/diff_peak/Act_KO_3_combine_peak.count.table", col_names = F)


# merge table -----------------------------------
colnames(count_df.ff_rep1) <- c("chrom", "start", "end", "Act_FF_1")
colnames(count_df.ff_rep2) <- c("chrom", "start", "end", "Act_FF_2")
colnames(count_df.ff_rep3) <- c("chrom", "start", "end", "Act_FF_3")

colnames(count_df.ko_rep1) <- c("chrom", "start", "end", "Act_KO_1")
colnames(count_df.ko_rep2) <- c("chrom", "start", "end", "Act_KO_2")
colnames(count_df.ko_rep3) <- c("chrom", "start", "end", "Act_KO_3")

count_df.merge <- count_df.ff_rep1

count_df.merge <- bind_cols(
  count_df.merge,
  dplyr::select(count_df.ff_rep2, Act_FF_2)
)

count_df.merge <- bind_cols(
  count_df.merge,
  dplyr::select(count_df.ff_rep3, Act_FF_3)
)

count_df.merge <- bind_cols(
  count_df.merge,
  dplyr::select(count_df.ko_rep1, Act_KO_1)
)


count_df.merge <- bind_cols(
  count_df.merge,
  dplyr::select(count_df.ko_rep2, Act_KO_2)
)

count_df.merge <- bind_cols(
  count_df.merge,
  dplyr::select(count_df.ko_rep3, Act_KO_3)
)


# save count table 
write_tsv(count_df.merge, file = "./result/downstream/diff_peak/merge_count_table.tsv")


# 测序深度在total_count.sh脚本中统计过，直接读取结果 -----------------------

# calculate total count 
case_info.vec = c(
  "Act_FF_1",
  "Act_FF_2",
  "Act_FF_3",
  "Act_KO_1",
  "Act_KO_2",
  "Act_KO_3"
)

merge_df <- NULL

for(case_info in case_info.vec){
  print(case_info)
  
  total_df_filename <- sprintf("./result/count_table/%s.total.count_table", case_info)  
  
  total_df <- read_tsv(file = total_df_filename, col_names = F)  
  colnames(total_df) <- c("chrom", "length", "map_read_count", "unmap_read_count")
  
  total_df$case_info = case_info
  
  
  merge_df <- bind_rows(
    merge_df,
    total_df
  )
  
}


# count total
total_count.vec = filter(
  merge_df,
  chrom != "chrY",
  chrom != "chrM"
) %>% group_by(
  case_info
) %>% summarise(
  total_count = sum(map_read_count)
) %>% pull(
  total_count
) 

# 计算得到的 total_count.vec 是从 merge_df 里按 case_info 分组汇总出来的 总 reads 数
total_count.vec

colnames(count_df.merge)

sum(merge_df$map_read_count[merge_df$case_info == "Act_FF_1"])
sum(merge_df$map_read_count[merge_df$case_info == "Act_FF_2"])
sum(merge_df$map_read_count[merge_df$case_info == "Act_FF_3"])
sum(merge_df$map_read_count[merge_df$case_info == "Act_KO_1"])
sum(merge_df$map_read_count[merge_df$case_info == "Act_KO_2"])
sum(merge_df$map_read_count[merge_df$case_info == "Act_KO_3"])

# 有时候顺序会反，可能跟字母相关！需要注意顺序
filter(
  merge_df,
  chrom != "chrY",
  chrom != "chrM"
) %>% group_by(
  case_info
) %>% summarise(
  total_count = sum(map_read_count)
)
#total_count.vec.order <- total_count.vec[c(3,4,1,2)]



# 归一化 (CPM) ----------------------------------------------------
count_df.merge.mat <- as.matrix(dplyr::select(count_df.merge, Act_FF_1:Act_KO_3))
# count_df.merge.mat 的列顺序往往是固定的样本顺序

count_df.merge.mat.norm = count_df.merge.mat / total_count.vec * 1e6

colnames(count_df.merge.mat.norm)

count_df.merge.mat.norm.df <- as.data.frame(
  count_df.merge.mat.norm
)


# 评估Repeat之间相关性 -------------------------------------

# rep1 vs rep2 
# make correlation plot
library(ggpmisc)
library(ggplot2)
p12 <- ggplot(count_df.merge.mat.norm.df, aes(x=log2(Act_FF_1+1), y = log2(Act_FF_2+1)) ) +
  geom_point(color="orange") +
  scale_fill_continuous() +
  theme_bw() + 
  xlim(0,10) + ylim(0,10) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),) + 
  geom_abline(intercept = 0,slope = 1,linetype="dashed",color="#333333") + 
  labs(title = "ChIP-seq rep1 vs rep2",x="Log2(CPM.Rep1 + 1)",y="Log2(CPM.Rep2 + 1)") + 
  # geom_smooth(method = "lm", formula = y~x) + 
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) #添加回归方程和调整R方 


p13 <- ggplot(count_df.merge.mat.norm.df, aes(x=log2(Act_FF_1+1), y = log2(Act_FF_3+1)) ) +
  geom_point(color="orange") +
  scale_fill_continuous() +
  theme_bw() + 
  xlim(0,10) + ylim(0,10) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),) + 
  geom_abline(intercept = 0,slope = 1,linetype="dashed",color="#333333") + 
  labs(title = "ChIP-seq rep1 vs rep2",x="Log2(CPM.Rep1 + 1)",y="Log2(CPM.Rep2 + 1)") + 
  # geom_smooth(method = "lm", formula = y~x) + 
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) #添加回归方程和调整R方 

p23 <- ggplot(count_df.merge.mat.norm.df, aes(x=log2(Act_FF_2+1), y = log2(Act_FF_3+1)) ) +
  geom_point(color="orange") +
  scale_fill_continuous() +
  theme_bw() + 
  xlim(0,10) + ylim(0,10) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),) + 
  geom_abline(intercept = 0,slope = 1,linetype="dashed",color="#333333") + 
  labs(title = "ChIP-seq rep1 vs rep2",x="Log2(CPM.Rep1 + 1)",y="Log2(CPM.Rep2 + 1)") + 
  # geom_smooth(method = "lm", formula = y~x) + 
  stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T) #添加回归方程和调整R方 


p12 + p13 + p23


# 多样本Peak CPM热图相关性 ----------------------
library(pheatmap)

mat <- count_df.merge.mat.norm.df %>%
  dplyr::select(Act_FF_1, Act_FF_2, Act_FF_3) %>%
  dplyr::mutate(across(everything(), ~ log2(.x + 1)))

cor_mat <- cor(mat, method="pearson")

pheatmap(cor_mat,
         display_numbers=TRUE,
         cluster_rows=FALSE, cluster_cols=FALSE,
         main="Replicate Correlation (Pearson)")




# foldchange 计算 (base CPM)-------------------------------------------
count_df.merge.mat.norm.df.AddFC <- mutate(
  count_df.merge.mat.norm.df,
  ff_mean = (Act_FF_1 + Act_FF_2 + Act_FF_3) / 3,
  ko_mean = (Act_KO_1 + Act_KO_2 + Act_KO_3) / 3,
  log2FC = log2(ko_mean / ff_mean)
)

hist(count_df.merge.mat.norm.df.AddFC$log2FC)

peak_id <- count_df.merge %>%
  transmute(peak = paste0(chrom, ":", start, "-", end)) %>%
  pull(peak)

length(peak_id)
nrow(count_df.merge.mat.norm.df.AddFC)

rownames(count_df.merge.mat.norm.df) <- peak_id

# DESeq2 (base count)---------------------------------------------
library(DESeq2)

colnames(count_df.merge.mat)

sample_df <- data.frame(
  condition = c(rep("ff",3), rep("ko",3)),
  cell_line = "293T"
)

des.obj <- DESeqDataSetFromMatrix(countData = count_df.merge.mat, colData = sample_df, design = ~condition)

des.obj.res <- DESeq(des.obj)

res.df <- as.data.frame(results(des.obj.res))

head(res.df)
nrow(res.df)

rownames(res.df) <- peak_id



# Diffbind ----------------------------------------------------

library(DiffBind)

# 第一步是读取一组峰集及其相关元数据。这可以通过样本表完成。读取峰集后，合并函数会找到所有重叠的峰，并导出覆盖所有提供的峰的唯一基因组区间集（即实验的共识峰集）。
# 如果一个区域出现在两个以上的样本中，则将其纳入共识集。该共识集代表了用于后续分析的候选结合位点的总体集合。

# R: 生成 DiffBind samplesheet（假设在 project root 运行）
bam_files <- list.files("./result/bam", pattern="sortp_genome_rmchrM_rmdup.bam$", full.names = TRUE)
peak_files <- list.files("./result/peak_calling", pattern="macs2_peaks.rmBlacklist.narrowPeak$", full.names = TRUE)

# 从文件名提取 sample name（去掉后缀）
bam_samples <- sub("\\.sortp_genome_rmchrM_rmdup\\.bam$", "", basename(bam_files))
peak_samples <- sub("\\.macs2_peaks.rmBlacklist.narrowPeak$", "", basename(peak_files))

# 简单匹配：以名字相等为配对（如 Act_FF_1）
common <- intersect(bam_samples, peak_samples)
if(length(common)==0) stop("Can't match BAM and Peaks by name — check filenames")

# 构建 dataframe（这里假设样本名里包含 Condition，可按需改）
df <- data.frame(
  SampleID = common,
  Tissue = NA,
  Factor = "ATAC",
  Condition = ifelse(grepl("KO", common), "KO", "Act"),
  Replicate = as.integer(gsub(".*_([0-9]+)$", "\\1", common)),  # 取末尾数字为 replicate（如 _1,_2）
  bamReads = file.path("./result/bam", paste0(common, ".sortp_genome_rmchrM_rmdup.bam")),
  bamControl = NA,
  Peaks = file.path("./result/peak_calling", paste0(common, ".macs2_peaks.rmBlacklist.narrowPeak")),
  PeakCaller = "narrow",
  stringsAsFactors = FALSE
)

# 若某些 sample 没找到 peak 文件，可在此打印警告并手动检查
missing_peaks <- setdiff(common, peak_samples)
if(length(missing_peaks)>0) warning("missing peaks for: ", paste(missing_peaks, collapse=", "))

# 写出 CSV
write.csv(df, "./result/downstream/diff_peak/samplesheet_diffbind.csv", row.names = FALSE, quote = FALSE)


samples <- read.csv('./result/downstream/diff_peak/samplesheet_diffbind.csv')
dbObj <- dba(sampleSheet=samples)

dbObj
# 6 Samples, 21630 sites in matrix (29170 total):
# 将所有样本的peak区域进行合并，得到29170个合并后的peak区域
# 在合并后的peak区域进行挑选，只保留这些区域中与至少两个（该数字可以通过 minOverlap 参数调整，默认为2）样本存在重叠的区域

# $allcalled 记录了每个样本在3795区域中是否有重叠的信息
sum(apply(dbObj$allcalled, 1, sum)>=2)

# Computing summits…：读入了所有样本（包括实验组和对照组）的bam文件，然后根据peaks区域内read的分布，找到每个样本peak的峰值（summit）所在位置。
# Re-centering peaks…：如果一个peak区域内的summit坐标，是根据所有包含样本summit坐标的均值（中位数）所决定的。
# 根据峰值，将peaks区域的长度统一，即以峰值为中心，向两侧拓展200bp范围（默认为200，可以通过 summits 参数进行修改），即此时的peaks长度统一为401bp
# 接着计算在拓展区域内，每个样本所包含的reads数目
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dbObj

##进行标准化（DBA_EDGER或DBA_DESEQ2两种方法），也可以不做
# dbObj <- dba.normalize(dbObj, method=DBA_EDGER, library=DBA_LIBSIZE_FULL, normalize=DBA_NORM_LIB)


# 样本在差异peaks区域内的read矩阵
dba.plotPCA(dbObj,  attributes=DBA_CONDITION, label=DBA_ID)

# 样本在peaks区域内的read矩阵
plot(dbObj) 
# 这里的数据不是原始的raw数据，而是经过Normalization处理后的矩阵，Normalization方法是简单的CPM等仅处理测序深度的方式



# 差异分析 ------------------------

# dba.contrast 默认要求每组样本起码有3个生物学重复
dbObj <- dba.contrast(dbObj, reorderMeta=list(Condition="Act"), minMembers = 2)

dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)

dba.show(dbObj, bContrasts=T)
# 对于差异结果的定义，默认是 "FDR<0.05 & |logFC|>0"，这个标准可以通过参数进行修改


plot(dbObj, contrast=1)

dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_CONDITION, label=DBA_ID)
# contrast 用于指定差异分析的结果，而数字1来自 dba.show(dbObj, bContrast=T) 结果，因为第1行代表了这次分析的结果，因此参数指定为1



dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dba.plotVenn(dbObj, contrast=1, bGain=T, bLoss=T, bAll = F)
dba.plotVenn(dbObj, contrast=1, bDB=T, bNotDB=F, bGain=T, bLoss=T, bAll=T)
# bDB：是否纳入差异结合位点，根据logFC以及FDR定义所得。
# bNotDB：是否纳入非差异结合位点，根据logFC以及FDR定义所得。
# bGain：是否纳入结合增强位点，仅根据logFC>0挑选
# bLoss：是否纳入结合减弱位点，仅根据logFC<0挑选
# bAll：是否纳入所有位点

dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)



dba.plotVolcano(dbObj,method=DBA_EDGER)


# 展示不同组别之间，在所有差异结合区域、上调区域以及下调区域中的结合信号强度分布
pvals <- dba.plotBox(dbObj)


dba.plotHeatmap(dbObj, correlations=F)


res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# Write to file
saveRDS(dbObj, file = './result/downstream/diff_peak/diffbind.rds')
out <- as.data.frame(res_deseq)
write.table(out, file="./result/downstream/diff_peak/diffbind_deseq2.txt", sep="\t", quote=F, row.names=F)











