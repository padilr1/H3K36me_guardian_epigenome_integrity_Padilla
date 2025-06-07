# plot genome browser tracks using pyGenomeTracks

## displaying regions of broad H3K9me3 loss in mMSC
```bash
pyGenomeTracks --tracks mMSC_PA_TKO_H3K9me3.ini --region chr18:14,354,887-53,302,122 --dpi 600 -o mMSC_PA_TKO_H3K9me3.pdf --height 3 --plotWidth 8 --trackLabelFraction 0 --fontSize 3.5
```

## displaying redistribution of H3K9me3 in mMSC
```bash
pyGenomeTracks --tracks mMSC_PA_TKO_H3K9me3_SUV39H1_HP1_ATACseq.ini --region chr5:57,999,380-73,643,774 --dpi 600 -o mMSC_PA_TKO_H3K9me3_SUV39H1_HP1_ATACseq.pdf --height 6.5 --plotWidth 8 --trackLabelFraction 0 --fontSize 3.5
```

## displaying regions of broad H3K9me3 loss in Cal27
```bash
pyGenomeTracks --tracks Cal27_PA_H3K36MOE_K9me3.ini --region chr10:73,863,533-123,034,084 --dpi 600 -o Cal27_PA_H3K36MOE_K9me3.pdf --height 4.5 --plotWidth 8 --trackLabelFraction 0.00000001 --trackLabelHAlign center --fontSize 3.5
```

# ChIPbinner analysis used to generate 100kb genome-wide scatterplots

## mMSC
```r
ChIPbinner::filter_low_counts(out_dir = "outdir",genome_assembly = "mm10",sample_bedfiles = c("mMSC.TKO.H3K9me3.100kb.bed","mMSC.PA.H3K9me3.100kb.bed"),cutoff = 100)

norm_bw(out_dir = "outdir",genome_assembly = "mm10",use_input = TRUE,depth_norm = FALSE,immunoprecipitated_binned_file = "mMSC.PA.H3K9me3.100kb_filt.bed",input_binned_file = "mMSC.PA_input.H3K9me3.100kb_filt.bed",pseudocount = 1e-1,raw_count_cutoff = 0)

norm_bw(out_dir = "outdir",genome_assembly = "mm10",use_input = TRUE,depth_norm = FALSE,immunoprecipitated_binned_file = "mMSC.TKO.H3K9me3.100kb_filt.bed",input_binned_file = "mMSC.TKO_input.H3K9me3.100kb_filt.bed",pseudocount = 1e-1,raw_count_cutoff = 0)

pre_clust(out_dir = "outdir",treated_samp_norm_bw = "mMSC.TKO.H3K9me3.100kb_filt.bw",wildtype_samp_norm_bw = "mMSC.PA.H3K9me3.100kb_filt.bw",output_filename = "mMSC_H3K9me3_WT_TKO_100kb")

density_based_scatterplot(out_dir = "outdir",genome_assembly = "mm10",treated_samp_norm_bw = "mMSC.TKO.H3K9me3.100kb.norm.bw",wildtype_samp_norm_bw = "mMSC.PA.H3K9me3.100kb.norm.bw",are_R_objects = FALSE,cell_line = "mMSC",histone_mark = "H3K9me3",annotated_clusters = "cons.mMSC.PA.TKO.H3K9me3.100kb.rda",number_of_clusters = 2,output_filename = "mMSC_H3K9me3_WT_TKO_100kb",title_of_plot = "H3K9me3 enrichment",plot_title_font_size = 8,pow = 0.7,legend_font_size = 5,min_x = -4.5,min_y = -4.5,max_x = 0,max_y = 0,hexbins = 50,show_scales = FALSE,xaxis_label = "mMSC PA",yaxis_label = "mMSC TKO",axis_title_font_size = 7,height_of_figure = 4.5,width_of_figure = 11,include_additional_density_plot = TRUE,filter_extreme_bins = TRUE)
```

## Cal27
```r
ChIPbinner::filter_low_counts(out_dir = "outdir",genome_assembly = "hg38",sample_bedfiles = c("Cal27.K36MOE.H3K9me3.100kb.bed","Cal27.PA.H3K9me3.100kb.bed"),cutoff = 100)

norm_bw(out_dir = "outdir",genome_assembly = "hg38",use_input = TRUE,depth_norm = FALSE,immunoprecipitated_binned_file = "Cal27.PA.H3K9me3.100kb_filt.bed",input_binned_file = "Cal27.PA_input.H3K9me3.100kb_filt.bed",pseudocount = 1e-1,raw_count_cutoff = 0)

norm_bw(out_dir = "outdir",genome_assembly = "hg38",use_input = TRUE,depth_norm = FALSE,immunoprecipitated_binned_file = "Cal27.K36MOE.H3K9me3.100kb_filt.bed",input_binned_file = "Cal27.K36MOE_input.H3K9me3.100kb_filt.bed",pseudocount = 1e-1,raw_count_cutoff = 0)

pre_clust(out_dir = "outdir",treated_samp_norm_bw = "Cal27.K36MOE.H3K9me3.100kb_filt.bw",wildtype_samp_norm_bw = "Cal27.PA.H3K9me3.100kb_filt.bw",output_filename = "Cal27_H3K9me3_WT_H3K36MOE_100kb")

density_based_scatterplot(out_dir = "outdir",genome_assembly = "hg38",treated_samp_norm_bw = "Cal27.H3K36M_OE.H3K9me3.100kb.norm.bw",wildtype_samp_norm_bw = "Cal27.PA.H3K9me3.100kb.norm.bw",are_R_objects = FALSE,cell_line = "Cal27",histone_mark = "H3K9me3",annotated_clusters = "cons.Cal27.PA.H3K36M_OE.H3K9me3.100kb.rda",number_of_clusters = 2,output_filename = "Cal27_H3K9me3_WT_H3K36MOE_100kb",title_of_plot = "H3K9me3 enrichment",plot_title_font_size = 8,pow = 0.7,legend_font_size = 5,min_x = -4.5,min_y = -4.5,max_x = 0,max_y = 0,hexbins = 50,show_scales = FALSE,xaxis_label = "Cal27 PA",yaxis_label = "Cal27 H3K36M-OE",axis_title_font_size = 7,height_of_figure = 4.5,width_of_figure = 11,include_additional_density_plot = TRUE,filter_extreme_bins = TRUE)
```

# merging 100kb bins
```bash
bedtools merge -i rawClustA.bed -d 1000000 > clusA.100kb.merged1mb.bed
bedtools merge -i rawClustB.bed -d 1000000 > clusB.100kb.merged1mb.bed

bedtools subtract -a chrom.bed -b clusB.100kb.merged1mb.bed > purported_clusA.bed
bedtools intersect -a purported_clusA.bed -b clusA.100kb.merged1mb.bed -wa > mMSC_H3K9me3_clusA.bed
mv clusB.100kb.merged1mb.bed mMSC_H3K9me3_clusB.bed
```

# generating heatmaps at clusA and clusB regions - example using mMSC
```bash
ref="mMSC_H3K9me3_clusA.bed mMSC_H3K9me3_clusB.bed"
ref_lab="clusA clusB" # for each reference, indicate a label, separated by a space
# signal samples
samp="mMSC_PA_merged.H3K9me3.ms_cpm.bw mMSC_TKO_merged.H3K9me3.ms_cpm.bw"
# sample labels
samp_lab="PA TKO"
# blacklist
blacklist_file="mm10_blacklist.bed"
computeMatrix reference-point -R ${ref} \
-S ${samp} \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50000 -a 100000 -b 100000 \
--missingDataAsZero --skipZeros \
-p 16 --samplesLabel ${samp_lab} --verbose \
-o ${mat}

plotHeatmap -m ${mat} \
-o ${plot} \
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


# aggregate plot of SUV39H1-Flag within genes
```bash
computeMatrix scale-regions -R silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed \
-S PA_SUV39H1_FLAG_inputNorm.bw TKO_SUV39H1_FLAG_inputNorm.bw \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel PA TKO --verbose \
-o ${mat}
# plot
plotProfile -m ${mat} -o ${outputGraphic}.aggregate_profile.pdf --dpi 600 --colors blue darkkhaki --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight 3.5 --plotWidth 4.2 --perGroup -y "Input-norm signal" --startLabel "TSS" --endLabel "TES" --legendLocation "best" --labelRotation 90 --plotFileFormat pdf
```


# plot heatmaps of SUV39H1-OE and SUV39H1-KD

## PA
```bash
ref="mMSC_H3K9me3_clusA.bed mMSC_H3K9me3_clusB.bed"
samp="mMSC_PA_SUV39H1_OE_MSnorm.bw mMSC_PA_SUV39H1_KD_MSnorm.bw"
# sample labels
samp_lab="PA_SUV39H1_OE PA_SUV39H1_KD"
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

## TKO
```bash
ref="mMSC_H3K9me3_clusA.bed mMSC_H3K9me3_clusB.bed"
samp="mMSC_TKO_SUV39H1_OE_MSnorm.bw mMSC_TKO_SUV39H1_KD_MSnorm.bw"
# sample labels
samp_lab="TKO_SUV39H1_OE TKO_SUV39H1_KD"
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

# percentage of H3K9me3 within clusters
```bash
plotEnrichment -b PA_H3K9me3.bam PA_SUV39H1_KD_H3K9me3.bam TKO_H3K9me3.bam TKO_SUV39H1_KD_H3K9me3.bam \
--BED clusA.bed clusB.bed \
--regionLabels "clusA" "clusB" \
-o ${output}.png \
--blackListFileName mm10_blacklist.bed \
--centerReads \
-p 12 \
--ignoreDuplicates \
-e 200 \
--minMappingQuality 5 \
--outRawCounts ${output}.txt
```

```r
raw <- fread("output.txt") %>% mutate(samp = gsub(".*/","",.$file)) %>% mutate(cond = gsub("_H3K9me3","",.$samp))

df <- raw %>% dplyr::select(c("cond","samp","featureType","percent")) 
# df$cond <- factor(x = df$cond,levels = c("PA-SUV39H1-WT","PA-SUV39H1-KD","TKO-SUV39H1-WT","TKO-SUV39H1-KD"))
df$cond <- factor(x = df$cond,levels = c("TKO-SUV39H1-KD","TKO","PA-SUV39H1-KD","PA"))
df$featureType <- factor(x=df$featureType,levels = c("clusA","clusB"))

pointSize=0.5

x = "cond"
y = "percent"

d <- df %>%
    dplyr::rename(x = !!x, y = !!y) %>%
  group_by(x,featureType) 

stats <- d %>%
    summarize(mean_y = mean(y),
              sd_y = sd(y))

geomBarWidth = 0.2
colors = c("darkorange","cornflowerblue")
fontsize = 5
ttl=""
ylab="% H3K9me3 signal"
  
ggplot(stats, aes(x=x, y=mean_y,fill=featureType)) + 
    geom_bar(stat = "identity", position = position_dodge(),show.legend = TRUE)+
    labs(x="",y=ylab,title=ttl)+ 
    scale_fill_manual(values=colors) +
    coord_flip()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=fontsize+1,family="Helvetica",colour = "black"),
          axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
          axis.text.x=element_text(size=fontsize,family="Helvetica",colour = "black"),
          plot.title = element_text(color = "black",size=fontsize+1,family="Helvetica",hjust=0.5,vjust=0.5),
          strip.background =element_rect(fill="white"),
          strip.text = element_text(
            size = fontsize + 1, color = "black",family = "Helvetica"),
          # legend
          legend.text=element_text(size=fontsize,family = "Helvetica",color = "black"),
          legend.title=element_blank(),
          legend.background = element_rect(fill="white",size = 1),
          legend.key=element_rect(fill="white",size = 0.1),
          legend.margin = unit(c(0,0,0,0),"mm"),
          legend.spacing = unit(c(0,0,0,0),"mm"),
          legend.position = "top",
          legend.justification = "right",
          legend.box.margin = margin(-1,-1,-1,-1),
          legend.key.size = unit(3,"mm"),
          plot.margin = unit(c(0.8,0.8,0.8,0.8), "mm"),
          panel.spacing = unit(0.5,'cm'),
          panel.spacing.y = unit(0.5,'cm'),
          panel.spacing.x = unit(0.5,'cm'))
ggsave(filename = "percentage_H3K9me3_clusAclusB.pdf",path = "outdir",device = "pdf",width = 2.5,height = 1.8,dpi = 600)
```



# SUV39H1 OE and KD within genes; metagene body plots

## PA
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
samples="mMSC_PA_SUV39H1_OE_MSnorm.bw mMSC_PA_SUV39H1_KD_MSnorm.bw"
sampleLabels="PA_SUV39H1_OE PA_SUV39H1_KD"
# cm
computeMatrix scale-regions -R ${regions} \
-S ${samples} \
-bl /lustre06/project/6007495/padilr1/genomes/blacklist/mm10/mm10.cut_and_run_blacklist.merged.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${sampleLabels} --verbose \
-o ${mat}
# plot
plotProfile -m ${mat} -o ${outputGraphic}.aggregate_profile.pdf --dpi 600 --colors darkorange black --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight 3.5 --plotWidth 4.2 --perGroup -y MS-norm signal --startLabel "TSS" --endLabel "TES" --legendLocation "best" --labelRotation 90 --plotFileFormat pdf
```
## TKO
```bash
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
samples="mMSC_TKO_SUV39H1_OE_MSnorm.bw mMSC_TKO_SUV39H1_KD_MSnorm.bw"
sampleLabels="TKO_SUV39H1_OE TKO_SUV39H1_KD"
# cm
computeMatrix scale-regions -R ${regions} \
-S ${samples} \
-bl /lustre06/project/6007495/padilr1/genomes/blacklist/mm10/mm10.cut_and_run_blacklist.merged.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${sampleLabels} --verbose \
-o ${mat}
# plot
plotProfile -m ${mat} -o ${outputGraphic}.aggregate_profile.pdf --dpi 600 --colors darkorange black --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight 3.5 --plotWidth 4.2 --perGroup -y MS-norm signal --startLabel "TSS" --endLabel "TES" --legendLocation "best" --labelRotation 90 --plotFileFormat pdf
```

# aggregate plot centered on TEs
```r
computeMatrix scale-regions -R TEs_in_clustA.bed \
-S ${samples} \
--samplesLabel ${sampleLabels} \
-bl mm10_blacklist.bed \
-b 300 -a 300 -m 300 -bs 10 \
--missingDataAsZero --skipZeros \
-p 6 \
--verbose \
-o ${mat}
Height=6
Width=7.5
plotProfile -m ${mat} -o ${output_graphic}.aggregate_profile.pdf --dpi 600 --colors blue darkkhaki  --numPlotsPerRow 1 --regionsLabel "TEs in active regions" --plotTitle "" --plotHeight ${Height} --plotWidth ${Width} --perGroup -y "MS-norm signal" --startLabel "start" --endLabel "end" --legendLocation "best" --labelRotation 90 --plotFileFormat pdf
```

# Detroit in genome browser tracks
```bash
pyGenomeTracks --tracks Detroit_PA_H3K36MOE_H3K9me3.ini --region chr10:52,296,181-86,097,426 --dpi 600 -o Detroit_PA_H3K36MOE_H3K9me3.pdf --height 4.5 --plotWidth 8 --trackLabelFraction 0.00000001 --trackLabelHAlign center --fontSize 3.5
```

# Detroit ChIPbinner
```r
Detroit_WT_H3K9me3_ms = 25.5039000/100
Detroit_H3K36MOE_H3K9me3_ms = 18.3273000/100
H3K36MOE_scaling_factor = Detroit_H3K36MOE_H3K9me3_ms/Detroit_WT_H3K9me3_ms

norm_bw(out_dir = "outdir",genome_assembly = "hg38",use_input = TRUE,depth_norm = TRUE,immunoprecipitated_binned_file = "Detroit_OE-K36M_ChIP_H3K9me3.10kb_filt.bed",input_binned_file = "Detroit_OE-K36M_input.10kb_filt.bed",pseudocount = 1,raw_count_cutoff = 0,scaling_factor = Detroit_H3K36MOE_H3K9me3_ms)

norm_bw(out_dir = "outdir",genome_assembly = "hg38",use_input = TRUE,depth_norm = TRUE,immunoprecipitated_binned_file = "Detroit562_PA_H3k9me3.10kb_filt.bed",input_binned_file = "Detroit562_PA_input.10kb_filt.bed",pseudocount = 1,raw_count_cutoff = 0,scaling_factor = Detroit_WT_H3K9me3_ms)

pre_clust(out_dir = "outdir",treated_samp_norm_bw = "Detroit_OE-K36M_ChIP_H3K9me3.10kb_filt.bw",wildtype_samp_norm_bw = "Detroit562_PA_H3k9me3.10kb_filt.bw",output_filename = "Detroit562_K36MOE_WT_H3K9me3_10kb",are_R_objects = FALSE)

annotate_clust(number_of_clusters = 2,matrix_file = "Detroit562_K36MOE_WT_H3K9me3_10kb_mat.csv",pooled_bed_file = "Detroit562_K36MOE_WT_H3K9me3_10kb_pooled.bed",hdbscan_output_file = "clus.Detroit562_K36MOE_WT_H3K9me3_10kb.1000.1000.txt",output_filename = "Detroit_H3K9me3_WT_H3K36MOE_10kb",out_dir = "outdir")

density_based_scatterplot(out_dir = "outdir",genome_assembly = "hg38",treated_samp_norm_bw = "Detroit_OE-K36M_ChIP_H3K9me3.10kb_filt.bw",wildtype_samp_norm_bw = "Detroit562_P16_H3k9me3.10kb_filt.bw",are_R_objects = FALSE,cell_line = "Detroit",histone_mark = "H3K9me3",annotated_clusters = "Detroit_H3K9me3_WT_H3K36MOE_10kb.annotated_clusters.rda",number_of_clusters = 2,output_filename = "Detroit_H3K9me3_WT_H3K36MOE_10kb",title_of_plot = "H3K9me3 enrichment",plot_title_font_size = 8,pow = 0.9,legend_font_size = 6,min_x = -2.5,min_y = -2.5,max_x = 0,max_y = 0,hexbins = 150,show_scales = FALSE,xaxis_label = "WT",yaxis_label = "H3K36M-OE",axis_title_font_size = 7,height_of_figure = 5,width_of_figure = 13,include_additional_density_plot = TRUE,filter_extreme_bins = FALSE)
```

# Detroit heatmaps
```bash
samp="Detroit562_H3K9me3_PA_ms_norm.bw Detroit562_H3K9me3_H3K36MOE_ms_norm.bw"
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

# profilePlyr analysis comparing ChIP-Seq signals in clusters A and B between two conditions using H3K9me3 as an example
```r
proplyrObject <- import_deepToolsMat("mMSC_PA_TKO_H3K9me3.mat.gz")
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
  labs(x="",y="log2(H3K9me3 signal + 1)") +
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
ggsave(filename = "mMSC_H3K9me3_signal_at_clusters_A_and_B.pdf",path = "outdir",width = 3.6,height=3.4,units="cm",dpi = 600,device = "pdf")
```