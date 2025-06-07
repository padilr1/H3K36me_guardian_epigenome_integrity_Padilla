# genome browser track of broadening of H3K27me3
```bash
pyGenomeTracks --tracks mMSC_PA_TKO_H3K9me3_H3K27me3_HP1.ini --region chr5:57,999,380-73,643,774 --dpi 600 -o output.pdf --height 4.5 --plotWidth 8 --trackLabelFraction 0 --fontSize 3.5
```

# heatmaps 

## H3K9me3
```bash
ref="mMSC_H3K9me3_clusA.bed mMSC_H3K9me3_clusB.bed"
samp="mMSC_PA_H3K9me3_MSnorm.bw mMSC_TKO_H3K9me3_MSnorm.bw"
# sample labels
samp_lab="PA TKO"
# computeMatrix
computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50000 -a 100000 -b 100000 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
# heatmap
plotHeatmap -m ${mat} \
-o ${plot}.pdf \
--dpi 600 --perGroup \
--colorMap "coolwarm" \
-T "" \
--regionsLabel "" "" \
--refPointLabel "" \
--verbose \
--whatToShow "heatmap and colorbar" \
--heatmapHeight 0.1 --heatmapWidth 1 \
--plotFileFormat pdf
```
## HP1
```bash
ref="mMSC_H3K9me3_clusA.bed mMSC_H3K9me3_clusB.bed"
samp="mMSC_PA_HP1_inputNorm.bw mMSC_TKO_HP1_inputNorm.bw"
# sample labels
samp_lab="PA TKO"
# computeMatrix
computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50000 -a 100000 -b 100000 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
# heatmap
plotHeatmap -m ${mat} \
-o ${plot}.pdf \
--dpi 600 --perGroup \
--colorMap "coolwarm" \
-T "" \
--regionsLabel "" "" \
--refPointLabel "" \
--verbose \
--whatToShow "heatmap and colorbar" \
--heatmapHeight 0.1 --heatmapWidth 1 \
--plotFileFormat pdf
```

## H3K27me3
```bash
ref="mMSC_H3K9me3_clusA.bed mMSC_H3K9me3_clusB.bed"
samp="mMSC_PA_H3K27me3_msNorm.bw mMSC_TKO_H3K27me3_msNorm.bw"
# sample labels
samp_lab="PA TKO"
# computeMatrix
computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50000 -a 100000 -b 100000 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
# heatmap
plotHeatmap -m ${mat} \
-o ${plot}.pdf \
--dpi 600 --perGroup \
--colorMap "coolwarm" \
-T "" \
--regionsLabel "" "" \
--refPointLabel "" \
--verbose \
--whatToShow "heatmap and colorbar" \
--heatmapHeight 0.1 --heatmapWidth 1 \
--plotFileFormat pdf
```

# agg plot of protein-coding genes for HP1
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
samples="PA_HP1_inputNorm.bw TKO_HP1_inputNorm.bw"
sampleLabels="PA_HP1 TKO_HP1"
# cm
computeMatrix scale-regions -R ${regions} \
-S ${samples} \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
--missingDataAsZero --skipZeros \
-p 16 --samplesLabel ${sampleLabels} --verbose \
-o ${mat}
# plot
plotProfile -m ${mat} -o ${outputGraphic}.aggregate_profile.pdf --dpi 600 --colors blue darkkhaki --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight 3.5 --plotWidth 4.2 --perGroup -y "Input-norm signal" --startLabel "TSS" --endLabel "TES" --legendLocation "best" --labelRotation 90 --plotFileFormat pdf
```

# genome browser tracks 

## mMSC
```bash
pyGenomeTracks --tracks mMSC_PA.TKO.clusB.new_ATAC_peaks.ini --region chr5:54,299,317-57,668,983 --dpi 600 -o mMSC_newATACseqPeaks_clusBregions.pdf --height 6 --plotWidth 8 --trackLabelFraction 0.00000001 --trackLabelHAlign center --fontSize 4.5
```

## Cal27
```bash
pyGenomeTracks --tracks Cal27_newATACpeaks_inClusB.ini --region chr5:155,125,727-157,121,563 --dpi 600 -o Cal27_newATACseqPeaks_clusBregions.pdf --height 6 --plotWidth 8 --trackLabelFraction 0.00000001 --trackLabelHAlign center --fontSize 4.5
```

# generate RPKM bigWigs for RNA-seq data
```bash
bamCoverage -b RNAseq.bam \
-o RPKM.bw \
-bl blacklist.bed \
-p 6 \
--normalizeUsing RPKM
```

# compute aggregate plots centered on new ATAC-seq peaks for ATAC, H3K27ac, H3K4me1 and RNA-seq in mMSC
```bash
ref="mMSC_TKO_ATACpeaks_within_clusB.bed"
ref_lab="ATAC" 
samp="PA.bw TKO.bw"
# sample labels
samp_lab="PA TKO"
# blacklist
blacklist_file="mm10_blacklist.bed"
# computeMatrix
computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50 -a 2500 -b 2500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
min=0
max=5
plotHeatmap -m ${mat} -o ${plot}.pdf --dpi 600 --colorMap Reds --verbose -x "" -y "" --regionsLabel "ATAC" --heatmapHeight 0.5 --heatmapWidth 2.0 --refPointLabel "center" --plotTitle "" --perGroup --legendLocation none --plotFileFormat pdf
```

# compute aggregate plots centered on new ATAC-seq peaks for ATAC, H3K27ac, H3K4me1 and RNA-seq in Cal27
```bash
ref="Cal27_H3K36MOE_ATACpeaks_within_clusB.bed"
ref_lab="ATAC" 
samp="PA.bw H3K36MOE.bw"
# sample labels
samp_lab="PA H3K36M-OE"
# blacklist
blacklist_file="hg38_blacklist.bed"
# computeMatrix
computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50 -a 2500 -b 2500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
min=0
max=5
plotHeatmap -m ${mat} -o ${plot}.pdf --dpi 600 --colorMap Reds --verbose -x "" -y "" --regionsLabel "ATAC" --heatmapHeight 0.5 --heatmapWidth 2.0 --refPointLabel "center" --plotTitle "" --perGroup --legendLocation none --plotFileFormat pdf
```

# profilePlyr analysis comparing HP1 in cluster A and B regions
```r
proplyrObject <- import_deepToolsMat("mMSC_PA_TKO_HP1.mat.gz")
proplyrObject_long <- profileplyr::summarize(proplyrObject, 
                                fun = rowMeans, 
                                output = "long") 
proplyrObject_long[1:3, ]
library(ggplot2)
ggplot(proplyrObject_long, aes(x = Sample, y = log(Signal))) + 
       geom_boxplot() +
  facet_wrap(~dpGroup)
proplyrObject_long <- proplyrObject_long %>% group_by(dpGroup, Sample)

clusters <- proplyrObject_long %>% `names<-`(c('cluster','range','sample','signal')) %>% group_by(cluster,sample)
stats <- compare_means(signal ~ sample,data=clusters,method="wilcox.test",group.by = "cluster",paired = TRUE) %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p",.$p.format))

fontsize=4
max_signal=max(log2(clusters$signal + 1))*1.05
ggplot(data = clusters,aes(x=sample,y=log2(signal + 1),fill=sample)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),size=0.1,linewidth=0.1) +
  geom_boxplot(width=0.2, color="black", alpha=0.5,outlier.size = 0.00001,outlier.alpha = 0.1,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("blue","darkkhaki")) +
  labs(x="",y="log2(HP1 signal + 1)") +
   geom_signif(comparisons =list(c("PA", "TKO")),textsize = 1.8,tip_length = 0.03,size=0.1,annotations=stats$p.signif[1],y_position = c(max_signal,max_signal))+
  facet_wrap(. ~ cluster) +
  geom_hline(yintercept = max_signal*1.5,colour="white") +
    theme(
      # panel
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis labels 
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black",hjust = 0.5,vjust=0.8,angle=45),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.ticks = element_line(size = 0.1),
      plot.title = element_text(hjust = 0.5,color = "white",size=fontsize,family="Helvetica"),
      # background
      strip.background =element_rect(fill="black"),
      strip.text = element_text(
        size = fontsize, color = "white",family = "Helvetica"),
      # legend
      legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
      legend.title=element_text(size=fontsize,family="Helvetica",color="black"),
      legend.background = element_rect(fill="white"),
      legend.key=element_rect(fill="white"),
      legend.position="right",
      legend.justification="right",
      legend.box.spacing = unit(-0.001, "cm"),
      # margins
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      panel.spacing = unit(1,'cm'),
      panel.spacing.y = unit(1,'cm'),
      panel.spacing.x = unit(0.1,'cm')) 
ggsave(filename = "mMSC_HP1_signal_at_clusters_A_and_B.pdf",path = "outdir",width = 3.6,height=3.4,units="cm",dpi = 600,device = "pdf")
```

# profilePlyr analysis comparing H3K27me3 signal in cluster A and B regions using mMSC as an example
```r
proplyrObject <- import_deepToolsMat("mMSC_PA_TKO_H3K27me3.mat.gz")
proplyrObject_long <- profileplyr::summarize(proplyrObject, 
                                fun = rowMeans, 
                                output = "long") 
proplyrObject_long[1:3, ]
library(ggplot2)
ggplot(proplyrObject_long, aes(x = Sample, y = log(Signal))) + 
       geom_boxplot() +
  facet_wrap(~dpGroup)
proplyrObject_long <- proplyrObject_long %>% group_by(dpGroup, Sample)

clusters <- proplyrObject_long %>% `names<-`(c('cluster','range','sample','signal')) %>% group_by(cluster,sample)
stats <- compare_means(signal ~ sample,data=clusters,method="wilcox.test",group.by = "cluster",paired = TRUE) %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p",.$p.format))

fontsize=4
max_signal=max(log2(clusters$signal + 1))*1.05
ggplot(data = clusters,aes(x=sample,y=log2(signal + 1),fill=sample)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),size=0.1,linewidth=0.1) +
  geom_boxplot(width=0.2, color="black", alpha=0.5,outlier.size = 0.00001,outlier.alpha = 0.1,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("blue","darkkhaki")) +
  labs(x="",y="log2(H3K27me3 signal + 1)") +
   geom_signif(comparisons =list(c("PA", "TKO")),textsize = 1.8,tip_length = 0.03,size=0.1,annotations=stats$p.signif[1],y_position = c(max_signal,max_signal))+
  facet_wrap(. ~ cluster) +
  geom_hline(yintercept = max_signal*1.5,colour="white") +
    theme(
      # panel
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis labels 
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black",hjust = 0.5,vjust=0.8,angle=45),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.ticks = element_line(size = 0.1),
      plot.title = element_text(hjust = 0.5,color = "white",size=fontsize,family="Helvetica"),
      # background
      strip.background =element_rect(fill="black"),
      strip.text = element_text(
        size = fontsize, color = "white",family = "Helvetica"),
      # legend
      legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
      legend.title=element_text(size=fontsize,family="Helvetica",color="black"),
      legend.background = element_rect(fill="white"),
      legend.key=element_rect(fill="white"),
      legend.position="right",
      legend.justification="right",
      legend.box.spacing = unit(-0.001, "cm"),
      # margins
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      panel.spacing = unit(1,'cm'),
      panel.spacing.y = unit(1,'cm'),
      panel.spacing.x = unit(0.1,'cm')) 
ggsave(filename = "mMSC_H3K27me3_signal_at_clusters_A_and_B.pdf",path = "outdir",width = 3.6,height=3.4,units="cm",dpi = 600,device = "pdf")
```

# HP1 signal aggregate plot centered on TEs
```r
computeMatrix scale-regions -R TEs_in_clustA.bed \
-S HP1_PA_inputNorm.bw HP1_TKO_inputNorm.bw \
--samplesLabel ${sampleLabels} \
-bl mm10_blacklist.bed \
-b 300 -a 300 -m 300 -bs 10 \
--missingDataAsZero --skipZeros \
-p 6 \
--verbose \
-o ${mat}
Height=6
Width=7.5
plotProfile -m ${mat} -o ${output_graphic}.aggregate_profile.pdf --dpi 600 --colors blue darkkhaki  --numPlotsPerRow 1 --regionsLabel "TEs in active regions" --plotTitle "" --plotHeight ${Height} --plotWidth ${Width} --perGroup -y "Input-norm signal" --startLabel "start" --endLabel "end" --legendLocation "best" --labelRotation 90 --plotFileFormat pdf
```

# heatmaps centered of H3K27me3 signal centered on clusA and clusB using Cal27 as an example
```bash
ref="Cal27_clusA.bed Cal27_clusB.bed"
ref_lab="clusA clusB" 
# signal samples
samp="Cal27_PA_H3K9me3_msNorm.bw Cal27_H3K36MOE_H3K9me3_msNorm.bw"
# sample labels
samp_lab="PA H3K36M-OE"
# computeMatrix
computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50000 -a 100000 -b 100000 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
# heatmap
plotHeatmap -m ${mat} \
-o ${plot}.pdf \
--dpi 600 --perGroup \
--colorMap "coolwarm" \
-T "" \
--regionsLabel "" "" \
--refPointLabel "" \
--verbose \
--whatToShow "heatmap and colorbar" \
--heatmapHeight 0.1 --heatmapWidth 1 \
--plotFileFormat pdf
```

# Detroit genome browser tracks
```bash
pyGenomeTracks --tracks Detroit_newATACpeaks_inClusB.ini --region chr5:2,935,983-4,514,870 --dpi 600 -o Detroit_newATACseqPeaks_clusBregions.pdf --height 6 --plotWidth 8 --trackLabelFraction 0.00000001 --trackLabelHAlign center --fontSize 4.5
```

# Detroit heatmaps centered on de novo ATAC-seq peaks, using ATAC-seq signal as an example
```bash
ref="Detroit_denovo_H3K36MEOE_ATACpeaksWithinClusB.bed"
# signal samples
samp="Detroit_PA_ATAC_merged_cpmNorm.bw Detroit_H3K36MOE_ATAC_merged_cpmNorm.bw"
# sample labels
samp_lab="PA H3K36M-OE"
# blacklist
blacklist_file="hg38_blacklist.bed"

computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50 -a 2500 -b 2500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
min=0
max=5
plotHeatmap -m ${mat} -o ${plot}_heatmap.pdf --dpi 600 --colorMap Reds --verbose -x "" -y "" --regionsLabel ${histoneMark} --heatmapHeight 0.5 --heatmapWidth 2 --refPointLabel "center" --plotTitle "" --perGroup --legendLocation none --plotFileFormat pdf --samplesLabel "" ""
```

# agg plot centered of DNA methylation on de novo ATACseq peaks
```bash
ref="mMSC_denovo_TKO_ATACpeaks_within_clusB.bed"
samp="WT_DNAme.bw TKO_DNAme.bw"
# sample labels
samp_lab="PA TKO"
# blacklist
blacklist_file="mm10_blacklist.bed"
# computeMatrix
computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50 -a 2500 -b 2500 \
-p 6 --samplesLabel ${samp_lab} --verbose \
-o ${mat}

plotProfile -m ${mat} -o ${plot}.pdf --dpi 600 --plotHeight ${plotHeight} --plotWidth ${plotWidth} --colors blue darkkhaki --refPointLabel "center" --yAxisLabel "% CpG methylation" --plotFileFormat "pdf" --perGroup --plotTitle "DNA methylation" --regionsLabel " "
```