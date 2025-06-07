# plotting genome-browser tracks using 'trackplot' 
```r
bigWigs = c("mMSC_QKO_merged_H3K36me2.bw",
            "mMSC_QuiKO_merged_H3K36me2.bw",
            "mMSC_QKO_merged_H3K27me3.bw",
            "mMSC_QuiKO_merged_H3K27me3")

bigWigs = read_coldata(bws = bigWigs, build = "mm10")

loci = "chr9:87,701,841-87,735,475"

t = track_extract(colData = bigWigs, loci = loci,build = "mm10",nthreads = 8,binsize = 50,padding = 10000)

track_plot(summary_list = t,col=c("darkorange","maroon","darkorange","maroon"),track_names = c("QKO-H3K36me2","QuiKO-H3K36me2","QKO-H3K27me3","QuiKO-H3K27me3"),draw_gene_track = TRUE,bw_track_height = 2,peaks_track_height = c(1),show_ideogram = FALSE,gene_track_height = 1,left_mar = 3,show_axis = TRUE,gene_fsize = 1.2,y_max = c(41,41,0.2,0.2),y_min = c(0,0,0,0))
```

# genes marked with residual H3K36me2 in QKO - log2FoldChanges
```r
QKO_H3K36me2_residualPeaks <- import.bed("mMSC_QKO_residual_H3K36me2_peaks.bed")
load("mm10_genes.RData")

QKO_H3K36me2_marked_genes <- IRanges::subsetByOverlaps(genes,QKO_H3K36me2_residualPeaks)

to_filt <- QKO_H3K36me2_marked_genes

load("mMSC_RNASeq_dds.RData")

dds$condition <- relevel(dds$condition,"QKO")

dds <- nbinomWaldTest(dds)
resultsNames <- resultsNames(dds) %>% as.data.frame()
resLFC <- lfcShrink(dds=dds,coef=6,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("gene_id") %>% na.omit()
# previous
resLFC_filt <- resLFC[resLFC$gene_id %in% to_filt$id,]
resLFC_filt$cond <- "QuiKO/QKO"
# QKO vs QuiKO
QKO_QuiKO_resLFC <- resLFC_filt
#### TKO vs QKO ####
dds$condition <- relevel(dds$condition,"TKO")
# model.matrix(~ condition
dds <- nbinomWaldTest(dds)
resultsNames <- resultsNames(dds) %>% as.data.frame()
resLFC <- lfcShrink(dds=dds,coef=2,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("gene_id") %>% na.omit()
# previous
resLFC_filt <- resLFC[resLFC$gene_id %in% to_filt$id,]
resLFC_filt$cond <- "QKO/TKO"
QKO_TKO_resLFC <- resLFC_filt
#### ASH1L-KO vs NSDunedit ####
dds$condition <- relevel(dds$condition,"NSD_unedit")
# model.matrix(~ condition
dds <- nbinomWaldTest(dds)
resultsNames <- resultsNames(dds) %>% as.data.frame()
resLFC <- lfcShrink(dds=dds,coef=4,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("gene_id") %>% na.omit()
# previous
resLFC_filt <- resLFC[resLFC$gene_id %in% to_filt$id,]
resLFC_filt$cond <- "ASH1L-KO/PA"
ASH1LKO_NSDunedit_resLFC <- resLFC_filt
#### combine dataframes ####
combined_resLFC <- rbind(ASH1LKO_NSDunedit_resLFC,QKO_QuiKO_resLFC,QKO_TKO_resLFC) %>% dplyr::select(c("cond","log2FoldChange"))
combined_resLFC$cond <- factor(x = combined_resLFC$cond,levels = c("QuiKO/QKO", "QKO/TKO", "ASH1L-KO/PA"))
# make ann text dataframe
ann_text <- data.frame(x=2,log2FoldChange =-0.5,lab = c("1","2","3"),
                       cond = factor(c("QuiKO/QKO", "QKO/TKO", "ASH1L-KO/PA"),levels = c("QuiKO/QKO", "QKO/TKO", "ASH1L-KO/PA")))

# plot
ggplot(data = combined_resLFC,aes(x=cond,y=log2FoldChange)) +
  geom_violinhalf(aes(fill=cond),show.legend = FALSE,position=position_nudge(x=0.1),size=0.2,linewidth=0.1) +
  geom_jitter(aes(color=cond),show.legend=FALSE,position=position_jitter(width=0.03, height=0.03),size=0.5,alpha=0.5) +
  geom_boxplot(width=0.2, color="black", alpha=0.3,outlier.size = 0.05,show.legend = FALSE,aes(color=cond,fill=cond),position=position_nudge(x=-0.3),linewidth=0.1) +
  scale_fill_manual(values=c("darksalmon","darkkhaki","limegreen")) +
  scale_colour_manual(values=c("darksalmon","darkkhaki","limegreen"))+
  labs(x="",y="log2FoldChange") +
  coord_flip() +
  facet_wrap(~cond,ncol=1,scales = "free_y")+
  geom_hline(yintercept =0,linetype="dashed",linewidth=0.1) +
  # text 
  geom_text(data=ann_text,label=c(nrow(QKO_QuiKO_resLFC[QKO_QuiKO_resLFC$log2FoldChange < 0,]),nrow(QKO_TKO_resLFC[QKO_TKO_resLFC$log2FoldChange < 0,]),nrow(ASH1LKO_NSDunedit_resLFC[ASH1LKO_NSDunedit_resLFC$log2FoldChange < 0,])),hjust=2,vjust=-2.8,size=1.9) +
  geom_text(data=ann_text,label=c(nrow(QKO_QuiKO_resLFC[QKO_QuiKO_resLFC$log2FoldChange > 0,]),nrow(QKO_TKO_resLFC[QKO_TKO_resLFC$log2FoldChange > 0,]),nrow(ASH1LKO_NSDunedit_resLFC[ASH1LKO_NSDunedit_resLFC$log2FoldChange > 0,])),hjust=-3,vjust=-2.8,size=1.9) +
  theme(
    # panel
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y.left  = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.line.x.bottom  = element_line(colour="black"),
    axis.line.x.top  = element_line(colour="black"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed",linewidth=0.1),
    # axis
    axis.title.y= element_text(size=7,family="Helvetica",colour = "black"),
    axis.title.x= element_text(size=7,family="Helvetica",colour = "black"),
    axis.text.x = element_text(size=7,family="Helvetica",colour = "black",hjust = 1,vjust=1),
    axis.text.y=element_text(size=7,family="Helvetica",colour = "black"),
    # labels
    strip.background =element_blank(),
    strip.text.x = element_blank(),
    strip.text = element_text(
      size = 7, color = "black",family = "Helvetica"),
    # legend
    legend.text=element_text(size=7,family = "Helvetica",color = "black"),
    legend.title=element_text(size=7,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    # margins
    plot.margin = unit(c(0.8, 0.8, 0.8, 0.8), "mm"),
    panel.spacing = unit(0.5,'cm'),
    panel.spacing.y = unit(0.5,'cm'),
    panel.spacing.x = unit(0.5,'cm'),
    legend.position="right",
    legend.justification="right",
    legend.box.spacing = unit(-0.001, "cm"))
```

# genes marked with residual H3K36me2 in QKO - FPKM
```r
QKO_H3K36me2_residualPeaks <- import.bed("mMSC_QKO_residual_H3K36me2_peaks.bed")
load("mm10_genes.RData")

QKO_H3K36me2_marked_genes <- IRanges::subsetByOverlaps(genes,QKO_H3K36me2_residualPeaks)

to_filt <- QKO_H3K36me2_marked_genes

load("mm10_genes.RData")
geneLengths <- fread("mm10_geneLengths.tsv")
mcols(dds)$basepairs <- geneLengths
norm.counts <- DESeq2::fpkm(dds,robust = TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("id") %>%
  mutate(ASH1L_KO_mean_exp=(rowMeans(.[grep("ASH1L_KO", names(.), value = TRUE)])))%>%
  mutate(QKO_mean_exp=(rowMeans(.[grep("QKO", names(.), value = TRUE)]))) %>%
  mutate(QuiKO_mean_exp=(rowMeans(.[grep("QuiKO", names(.), value = TRUE)]))) %>%
  mutate(TKO_mean_exp=(rowMeans(.[grep("TKO", names(.), value = TRUE)]))) %>%
  mutate(NSD3KO_mean_exp=(rowMeans(.[grep("NSD3_KO", names(.), value = TRUE)]))) %>%
  mutate(PA_mean_exp=(rowMeans(.[grep("NSD1_KO|sgNSD3_unedited", names(.), value = TRUE)]))) %>% dplyr::select(c("id","PA_mean_exp","ASH1L_KO_mean_exp","TKO_mean_exp","QKO_mean_exp","QuiKO_mean_exp"))

norm_counts_filt <- norm.counts[norm.counts$id %in% to_filt$id,]

norm_counts_filt_longer <- norm_counts_filt %>% pivot_longer(cols = c(2:6),names_to = "cond",values_to="FPKM")
norm_counts_filt_longer$cond <- factor(norm_counts_filt_longer$cond,levels = c("PA_mean_exp","ASH1L_KO_mean_exp","TKO_mean_exp","QKO_mean_exp","QuiKO_mean_exp"))

norm_counts_filt_longer_not_logged <- norm_counts_filt_longer

get_box_stats <- function(y, upper_limit = max(norm_counts_filt_longer$FPKM) * 1.25) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}
stats_QKO_QuiKO <- compare_means(FPKM ~ cond,data=norm_counts_filt_longer_not_logged,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p=",.$p.format)) %>% dplyr::slice(10)
stats_TKO_QKO <- compare_means(FPKM ~ cond,data=norm_counts_filt_longer_not_logged,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p=",.$p.format)) %>% dplyr::slice(8)
stats_PA_ASH1LKO <- compare_means(FPKM ~ cond,data=norm_counts_filt_longer_not_logged,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p=",.$p.format)) %>% dplyr::slice(1)
# ggplot2
ggplot(data = norm_counts_filt_longer,aes(x=cond,y=FPKM,fill=cond)) +
  geom_violin(show.legend = FALSE,linewidth=0.2) +
  geom_boxplot(width=0.2, color="black", alpha=0.1,outlier.size = 0.01,show.legend = FALSE,linewidth=0.1) +
  scale_x_discrete(labels=c("PA","ASH1L-KO","TKO","QKO","QuiKO")) +
  scale_fill_manual(values=c("royalblue1","limegreen","darkkhaki","darkorange","maroon")) +
  labs(x="",y="log2(FPKM+1)") +
  # signif annotations
  geom_signif(textsize = 2.5,tip_length = 0.01,xmin=stats_QKO_QuiKO$group1,xmax=stats_QKO_QuiKO$group2,annotations=stats_QKO_QuiKO$p.signif,y_position = max(norm_counts_filt_longer$FPKM) * 1.05,size=0.3) +
  geom_signif(textsize = 2.5,tip_length = 0.01,xmin=stats_TKO_QKO$group1,xmax=stats_TKO_QKO$group2,annotations=stats_TKO_QKO$p.signif,y_position = max(norm_counts_filt_longer$FPKM) * 1.03,size=0.3) +
  geom_signif(textsize = 2.5,tip_length = 0.01,xmin=stats_PA_ASH1LKO$group1,xmax=stats_PA_ASH1LKO$group2,annotations=stats_PA_ASH1LKO$p.signif,y_position = max(norm_counts_filt_longer$FPKM) * 1.01,size=0.3) +
  geom_hline(yintercept = 7.5,colour="white") +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y= element_text(size=7,family="Helvetica",colour = "black"),
    # angle = 45,hjust = 1,vjust=1
    axis.text.x = element_text(size=7,family="Helvetica",colour = "black"),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=7,family="Helvetica",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "blue",size=7,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
      size = 7, color = "black",family = "Helvetica"),
    legend.text=element_text(size=7,family = "Helvetica",color = "black"),
    legend.title=element_text(size=7,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
    panel.spacing = unit(0.5,'cm'),
    panel.spacing.y = unit(0.5,'cm'),
    panel.spacing.x = unit(0.5,'cm'),
    legend.position="right",
    legend.justification="right",
    legend.box.spacing = unit(-0.001, "cm"))
```

# correlation between H3K36me2 and gene expression changes
```r
QKO_H3K36me2_residualPeaks <- import.bed("mMSC_QKO_residual_H3K36me2_peaks.bed")
load("mm10_genes.RData")

QKO_H3K36me2_marked_genes <- IRanges::subsetByOverlaps(genes,QKO_H3K36me2_residualPeaks)

to_filt <- QKO_H3K36me2_marked_genes

load("mm10_genes.RData")
counts <- fread("mMSC_QKO_H3K36me2_epic2_peaks.counts")
counts$QKO_rep1 <- ((counts$QKO_rep1)/16947864)*1e6
counts$QKO_rep2 <- ((counts$QKO_rep2)/5778384)*1e6
counts$QKO_rep3 <- ((counts$QKO_rep3)/18895221)*1e6
# mutate to get row average 
counts <- counts %>% as.data.frame() %>% mutate(QKO_mean_peak_CPM=(rowMeans(.[grep("QKO", names(.), value = TRUE)])))
# write bed file
counts_gr <- counts[,c("Chr","Start","End","Geneid","QKO_mean_peak_CPM","Strand")] %>% makeGRangesFromDataFrame(.,keep.extra.columns = TRUE,seqnames.field = "Chr",start.field = "Start",end.field = "End",strand.field = "Strand")

QKO_H3K36me2_marked_genes <- counts_gr %>% IRanges::subsetByOverlaps(genes,counts_gr)

load("mMSC_RNASeq_dds.Rdata")
# run log2FoldChange shrinkage
dds$condition <- relevel(dds$condition,"QKO")
# model.matrix(~ condition
dds <- nbinomWaldTest(dds)
resultsNames <- resultsNames(dds) %>% as.data.frame()
resLFC <- lfcShrink(dds=dds,coef=6,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("gene_id") %>% na.omit()
# previous
resLFC_filt <- resLFC[resLFC$gene_id %in% to_filt$id,]
resLFC_filt$cond <- "QuiKO/QKO"
genes_df <- genes %>% as.data.frame()
genes_df$gene_id <- gsub("\\..*","",genes_df$gene_id)
QKO_genes <- resLFC_filt %>% left_join(genes_df,by="gene_id")
QKO_gene_bed <- QKO_genes[,c("seqnames","start","end","gene_name","log2FoldChange","strand")]
write_tsv(QKO_gene_bed,file="QKO_genes.bed",col_names = FALSE,quote = "none")
```

```bash
bedtools intersect -a QKO_genes.bed -b mMSC_QKO_residual_H3K36me2_peaks.bed -wao > QKO_intersect_peaks_genes.bed
```

```r
agg <- fread("QKO_intersect_peaks_genes.bed")

ggplot(data=agg,aes(x=V5,y=log2(V11))) +
  geom_point(size=0.7,alpha=0.8) +
  stat_cor(method="pearson",label.x = -1.5,label.y = 13,size=1.97) +
  labs(x="Gene Expression log2 (FC)\n QuiKO/QKO",y="QKO H3K36me2 log2(CPM)") +
  geom_smooth(method=lm,se=FALSE,fullrange=TRUE) +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=0.2),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y= element_text(size=6,family="Helvetica",colour = "black"),
    axis.title.x= element_text(size=6,family="Helvetica",colour = "black"),
    axis.text.x = element_text(size=6,family="Helvetica",colour = "black",hjust = 1,vjust=1),
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=6,family="Helvetica",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "blue",size=6,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
      size = 6, color = "black",family = "Helvetica"),
    legend.text=element_text(size=6,family = "Helvetica",color = "black"),
    legend.title=element_text(size=6,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "mm"),
    panel.spacing = unit(0.5,'cm'),
    panel.spacing.y = unit(0.5,'cm'),
    panel.spacing.x = unit(0.5,'cm'),
    legend.position="right",
    legend.justification="right",
    legend.box.spacing = unit(-0.001, "cm"))
```


# compute accessible enhancers
```r
QKO_ATAC <- fread("QKO_ATAC_IDR_peaks.bed") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
QuiKO_ATAC <- fread("QuiKO_ATAC_IDR_peaks.bed") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
# intersect 
ATAC_ol <- QKO_ATAC[overlapsAny(QKO_ATAC,QuiKO_ATAC)] %>% GenomicRanges::reduce()
write_tsv(ATAC_ol %>% as.data.frame() %>% dplyr::select(c(1:3)),file="QKO_QuiKO_intersect_ATAC_peaks.bed",col_names = FALSE,quote = "none")
# mm10 genes
mm10_genes <- fread("mm10_genes.bed")
promoters <- promoters(mm10_genes,upstream = 3000,downstream = 3000)
# find ATAC-seq not in promoter regions
accessible_enhancers <- ATAC_ol[!overlapsAny(ATAC_ol,promoters)] %>% GenomicRanges::reduce()
```

# compute accessible promoters
```r
QKO_ATAC <- fread("QKO_ATAC_IDR_peaks.bed") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
QuiKO_ATAC <- fread("QuiKO_ATAC_IDR_peaks.bed") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
# intersect 
ATAC_ol <- QKO_ATAC[overlapsAny(QKO_ATAC,QuiKO_ATAC)] %>% GenomicRanges::reduce()
write_tsv(ATAC_ol %>% as.data.frame() %>% dplyr::select(c(1:3)),file="QKO_QuiKO_intersect_ATAC_peaks.bed",col_names = FALSE,quote = "none")
# mm10 genes
mm10_genes <- fread("mm10_genes.bed")
promoters <- promoters(mm10_genes,upstream = 1500,downstream = 500)
accessible_promoters <- ATAC_ol[overlapsAny(ATAC_ol,promoters)] %>% GenomicRanges::reduce()
```

# plot PTM, ATAC and DNAme signal centered on accessible enhancers and promoters
```bash
accessibleRegions="QKO_QuiKO_accessible_enhancers.bed"
accessibleRegions="QKO_QuiKO_accessible_promoters.bed"
computeMatrix reference-point -R ${accessibleRegions} \
-S QKO.bw QuiKO.bw \
-bl mm10_blacklist.bed \
--referencePoint center \
--binSize 50 -a 2500 -b 2500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel QKO QuiKO --verbose \
-o out.mat.gz

plotProfile -m out.mat.gz -o aggregate_plot.pdf --dpi 600 --plotHeight 3 --plotWidth 3 --colors darkorange maroon --refPointLabel "center" --yAxisLabel "" --plotFileFormat "pdf" --perGroup --regionsLabel " "
```

# agg plot of PA and ASH1L-KO on promoters
```bash
computeMatrix reference-point -R promoters_marked_by_residual_H3K36me2_in_QKO.bed \
-S PA_H3K36me2_msNorm.bw ASH1LKO_H3K36me2_msNorm.bw \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50 -a 3000 -b 3000 \
-p 16 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
# plot
plotHeight=3.5
plotWidth=3.5
colour_list="blue limegreen"
plotProfile -m ${mat} -o ${plot}.pdf --dpi 600 --plotHeight ${plotHeight} --plotWidth ${plotWidth} --colors ${colour_list} --refPointLabel "TSS" --yAxisLabel "" --plotFileFormat "pdf" --perGroup --plotTitle "H3K36me2" --regionsLabel " "
```