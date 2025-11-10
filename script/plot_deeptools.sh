
#mkdir ../result/downstream/deeptools

# plot with TSS
# computeMatrix reference-point \
# -S \
# ../result/bigwig/Act_FF_1.RPKM.bw \
# ../result/bigwig/Act_KO_1.RPKM.bw \
# -R \
# ../ref_data/mm10/ataqc/mm10.tss.bed \
# --referencePoint center \
# --beforeRegionStartLength 2000 \
# --afterRegionStartLength 2000 \
# --binSize 10 \
# --skipZeros \
# --numberOfProcessors 20 \
# -o ../result/downstream/deeptools/tss_enh_region_Act_sample.profile.mat.gz \
# --samplesLabel Act_FF  Act_KO

# plot with TSS
plotHeatmap -m ../result/downstream/deeptools/tss_enh_region_Act_sample.profile.mat.gz \
    --colorMap Purples \
    -out ../result/downstream/deeptools/tss_enh_region_Act_sample.profile.heatmap.pdf



