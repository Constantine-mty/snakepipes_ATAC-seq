gunzip -c ENCFF493CCB.bed.gz > hg38.tss.bed
gunzip -c ENCFF356LFX.bed.gz > hg38.blacklist.bed


# 生成上下游各 2kb 的TSS扩展区域
awk 'BEGIN{OFS="\t"}{
    chr=$1;
    start=$2;
    end=$3;
    gene=$4;
    score=$5;
    strand=$6;

    if(strand=="+"){
        tss=start;
    } else if(strand=="-"){
        tss=end;
    }

    new_start = tss - 2000;
    new_end   = tss + 2000;
    if(new_start < 0) new_start = 0;

    print chr, new_start, new_end, gene, score, strand
}' hg38.tss.bed \
| sort -k1V -k2,2n > hg38.tss_extend2kb.sorted.bed


echo -e "blacklist_unzip\t${PWD}/hg38.blacklist.bed" >> ../hg38.tsv

echo -e "tss_extend2kb\t${PWD}/hg38.tss_extend2kb.sorted.bed" >> ../hg38.tsv

