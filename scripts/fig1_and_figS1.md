# plotting genome-browser tracks using pyGenomeTracks 

## top half
```bash
make_tracks_file --trackFiles mMSC_PA_merged_H3K36me2.bw mMSC_DKO_merged_H3K36me2.bw mMSC_TKO_merged_H3K36me2.bw mMSC_QKO_merged_H3K36me2.bw mm10.ncbiRefSeq.gtf -o ${track}.ini

pyGenomeTracks --tracks ${track}.ini --region chr5:53400510-53644953 --dpi 1200 -o ${track}.pdf --height 5.5 --plotWidth 9.5 --trackLabelFraction 0 --fontSize 5
```

## bottom half
```bash
make_tracks_file --trackFiles mMSC_TKO_merged_H3K36me2.bw mMSC_QKO_merged_H3K36me2.bw mMSC_TKO_merged_ATAC.bw mMSC_QKO_merged_ATAC.bw mMSC_TKO_merged_H3K27ac.bw mMSC_QKO_merged_H3K27ac.bw active_enhancers.bed strong_enhancers.bed mm10.ncbiRefSeq.gtf -o ${track}.ini

pyGenomeTracks --tracks ${track}.ini --region chr5:53,294,411-54,154,525 --dpi 1200 -o ${track}.pdf --height 5.5 --plotWidth 9.5 --trackLabelFraction 0 --fontSize 5
```

# compute accessible enhancers
```r
TKO_ATAC <- fread("TKO_ATAC_IDR_peaks.bed") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
QKO_ATAC <- fread("QKO_ATAC_IDR_peaks.bed") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
# intersect 
ATAC_ol <- TKO_ATAC[overlapsAny(TKO_ATAC,QKO_ATAC)] %>% GenomicRanges::reduce()
write_tsv(ATAC_ol %>% as.data.frame() %>% dplyr::select(c(1:3)),file="TKO_QKO_intersect_ATAC_peaks.bed",col_names = FALSE,quote = "none")
# mm10 genes
mm10_genes <- fread("mm10_genes.bed")
promoters <- promoters(mm10_genes,upstream = 3000,downstream = 3000)
# find ATAC-seq not in promoter regions
accessible_enhancers <- ATAC_ol[!overlapsAny(ATAC_ol,promoters)] %>% GenomicRanges::reduce()
```

# compute strong enhancers
```bash
ROSE_main.py -g mm10 -i TKO_QKO_accessible_enhancers.gff -r mMSC_TKO_H3K27ac_ChIP.bam -o outdir -t 2500 -c mMSC_TKO_input.bam
```

# compute accessible promoters
```r
TKO_ATAC <- fread("TKO_ATAC_IDR_peaks.bed") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
QKO_ATAC <- fread("QKO_ATAC_IDR_peaks.bed") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
# intersect 
ATAC_ol <- TKO_ATAC[overlapsAny(TKO_ATAC,QKO_ATAC)] %>% GenomicRanges::reduce()
write_tsv(ATAC_ol %>% as.data.frame() %>% dplyr::select(c(1:3)),file="TKO_QKO_intersect_ATAC_peaks.bed",col_names = FALSE,quote = "none")
# mm10 genes
mm10_genes <- fread("mm10_genes.bed")
promoters <- promoters(mm10_genes,upstream = 1500,downstream = 500)
accessible_promoters <- ATAC_ol[overlapsAny(ATAC_ol,promoters)] %>% GenomicRanges::reduce()
```

# plot PTM, ATAC and DNAme signal centered on accessible enhancers and promoters
```bash
accessibleRegions="TKO_QKO_accessible_enhancers.bed"
accessibleRegions="TKO_QKO_accessible_promoters.bed"
computeMatrix reference-point -R ${accessibleRegions} \
-S TKO.bw QKO.bw \
-bl mm10_blacklist.bed \
--referencePoint center \
--binSize 50 -a 2500 -b 2500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel TKO QKO --verbose \
-o out.mat.gz

plotProfile -m out.mat.gz -o aggregate_plot.pdf --dpi 600 --plotHeight 3 --plotWidth 3 --colors olive darkorange --refPointLabel "center" --yAxisLabel "" --plotFileFormat "pdf" --perGroup --regionsLabel " "
```

# plot PTM, ATAC and DNAme signal centered on strong enhancers
```bash
computeMatrix scale-regions -R strongEnhancers.bed \
-S TKO.bw QKO.bw \
-bl mm10_blacklist.bed \
--binSize 200 -a 5000 -b 5000 -m 10000 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel TKO QKO --verbose \
-o out.mat.gz

plotProfile -m out.mat.gz -o aggregatePlot.pdf --dpi 600 --plotHeight 3 --plotWidth 3 --colors olive darkorange --yAxisLabel "" --plotFileFormat "pdf" --perGroup --regionsLabel " " --startLabel "Start" --endLabel "End" --labelRotation 90
```

# compare genes within large H3K36me2 peaks and genes outside of large H3K36me2 peaks

```bash
# compute H3K36me2 peaks for each TKO replicate using epic2
epic2 --treatment mMSC_TKO_H3K36me2.bam \
--control mMSC_TKO_input.bam \
--genome mm10 \
-fdr 0.01 \
--output mMSC_TKO_H3K36me2_epic2.peaks
```

```bash
# run featureCounts on peaks
featureCounts -a mMSC_TKO_H3K36me2_epic2_peaks.saf" -F SAF -o mMSC_TKO_H3K36me2_peaks.counts -T 6 -p mMSC_H3K36me2_TKO_rep1.bam mMSC_H3K36me2_TKO_rep2.bam mMSC_H3K36me2_TKO_rep3.bam
```
## compute large TKO H3K36me2 peaks
```r
counts <- fread("mMSC_TKO_H3K36me2_peaks.counts")
# normalize raw counts by library size
counts$mMSC_H3K36me2_TKO_rep1 <- ((counts$mMSC_H3K36me2_TKO_rep1)/56480636)*1e6
counts$mMSC_H3K36me2_TKO_rep2 <- ((counts$mMSC_H3K36me2_TKO_rep2)/62150184)*1e6
counts$mMSC_H3K36me2_TKO_rep3 <- ((counts$mMSC_H3K36me2_TKO_rep3)/50021455)*1e6
# mutate to get row average 
counts <- counts %>% as.data.frame() %>% mutate(TKO_mean_peak_CPM=(rowMeans(.[grep("TKO", names(.), value = TRUE)]))) %>% dplyr::filter(TKO_mean_peak_CPM >= 100)
counts_gr <- counts[,c("Chr","Start","End","Geneid","TKO_mean_peak_CPM","Strand")] %>% makeGrangesFromDataframe
# write_tsv(counts_bed,file="~/Documents/10T_downstream/work/TKO_QKO_expression_K36me2_analyses/TKO_peaks.bed",col_names = FALSE,quote = "none")
# process peaks with at least 100
counts_filt <- counts_proc %>% dplyr::filter(TKO_mean_peak_CPM >= 100)
counts_gr <- counts_filt[,c("Chr","Start","End","Geneid","TKO_mean_peak_CPM","Strand")] %>% makeGRangesFromDataFrame(.,keep.extra.columns = TRUE,seqnames.field = "Chr",start.field = "Start",end.field = "End",strand.field = "Strand")
# import enhancers from FANTOM and Ensembl
mm10_ENSEMBL_enhancers <- import.bed("mm10EnsemblEnhancers.bed")
mm10_FANTOM_enhancers <- import.bed("mm10FantomEnhancers.bed")

mm10_annotated_enhancers <- c(mm10_ENSEMBL_enhancers,mm10_FANTOM_enhancers) %>% GenomicRanges::reduce() 

counts_gr <- counts_gr %>% IRanges::subsetByOverlaps(counts_gr,mm10_annotated_enhancers)

large_H3K36me2_enhancers_peaks <- counts_gr + 50000
```

## filter for genes with at least 1 FPKM 
```r
mMSC_TKO_expressedGenes <- mMSC_FPKM_counts %>% dplyr::filter(TKO_FPKM > 1)
filtered_mm10_protein_coding_genes <- dplyr::filter(Geneid %in% mMSC_TKO_expressedGenes$Geneid)
write_tsv(filtered_mm10_protein_coding_genes,file="filtered_mm10_protein_coding_genes.bed",col_names = FALSE,quote = "none")
```

## identify genes within large H3K36me2 enhancer peaks
```bash
bedtools intersect -a filtered_mm10_protein_coding_genes.bed -b mMSC_large_TKO_H3K36me2_enhancer_peaks.bed -wao > proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks.bed
```

## compute expression-matched genes and prepare matrix for plotting
```r
proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks <- fread("proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks.bed") 
colnames(proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks) <- c("chrom","TSS","TES","Geneid","log2FoldChange","strand","peak_chrom","start","end","peak_id","peak_count","strand","number_of_overlaps")
proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks <- proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks %>% dplyr::filter(peak_count != -1) %>% left_join(.,QKO_TKO_resLFC,by="Geneid") %>% dplyr::mutate(type="Genes within large H3K36me2 peaks") %>% dplyr::select(c("type","log2FoldChange"))

filt_FPKM_counts <- mMSC_FPKM_counts %>%
  mutate(TKO_mean_exp=(rowMeans(.[grep("TKO", names(.), value = TRUE)]))) %>% dplyr::select(c("id","TKO_mean_FPKM")) %>%
  dplyr::filter(TKO_mean_FPKM > 1) %>%
  mutate(bin=cut(TKO_mean_FPKM,breaks=seq(0,2050,by=3)))
expressionMatchedGenes <- filt_FPKM_counts %>% dplyr::filter(TKO_mean_FPKM >= 9 & TKO_mean_FPKM <= 12) %>% dplyr::filter(!Geneid %in% proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks$Geneid) %>% left_join(.,QKO_TKO_resLFC,by="Geneid") %>% dplyr::slice_sample(n=nrow(proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks),replace=FALSE) %>% dplyr::mutate(type="Genes outside of large H3K36me2 peaks") %>% dplyr::select(c("type","log2FoldChange"))

agg <- rbind(proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks,expressionMatchedGenes)
```

## plot
```r
size_of_text <- as.numeric(2)
stats <- compare_means(log2FoldChange ~ type,data=agg,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("****"))
# run ggplot
ggplot(data = agg,aes(x=type,y=log2FoldChange,fill=type)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1)) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.1) +
  geom_boxplot(width=0.1, color="black", alpha=0.1,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2)) +
  scale_fill_manual(values=c("burlywood","brown3")) +
  scale_x_discrete(labels=c("Genes within\n large H3K36me2\n peaks","Genes outside of\n large H3K36me2\n peaks")) +
  geom_hline(yintercept =0,linetype="dashed") +
  geom_hline(yintercept = 2,color="white")+
  labs(x="",y="Log2 (FC)\ngene expression QKO/TKO") +
  # stats
  geom_signif(textsize = 2,tip_length = 0.004,xmin=stats$group1,xmax=stats$group2,annotations=stats$pvalue,y_position = max(agg$log2FoldChange * 1.05)) +
  # count the number that go up
  annotate("text",y=0.7, x = 1.5, label = paste0(nrow(proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks[proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks$log2FoldChange > 0,])),size=size_of_text) +
  annotate("text",y=0.7, x = 2.5, label = paste0(nrow(expressionMatchedGenes[expressionMatchedGenes$log2FoldChange> 0,])),size=size_of_text) +
  # count those that go down
   annotate("text",y=-0.7, x = 1.5, label = paste0(nrow(proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks[proteinCodingGenes_within_TKOlargeH3K36me2EnhancerPeaks$log2FoldChange < 0,])),size=size_of_text) +
  annotate("text",y=-0.7, x = 2.5, label = paste0(nrow(expressionMatchedGenes[expressionMatchedGenes$log2FoldChange < 0,])),size=size_of_text) +
  coord_flip() +
  theme(
    # panel
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y.left  = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.line.x.bottom  = element_line(colour="black"),
    axis.line.x.top  = element_line(colour="black"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
    plot.title = element_text(hjust = 0.5,color = "black",size=6,family="Helvetica"),
    # axis
    axis.title.y= element_text(size=6,family="Helvetica",colour = "black"),
    axis.title.x= element_text(size=6,family="Helvetica",colour = "black"),
    axis.text.x = element_text(size=6,family="Helvetica",colour = "black",hjust = 1,vjust=1),
    axis.text.y=element_text(size=6,family="Helvetica",colour = "black"),
    # labels
    strip.background =element_blank(),
    strip.text.x = element_blank(),
    strip.text = element_text(
      size = 6, color = "black",family = "Helvetica"),
    # legend
    legend.text=element_text(size=6,family = "Helvetica",color = "black"),
    legend.title=element_text(size=6,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    # margins
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
    panel.spacing = unit(0.5,'cm'),
    panel.spacing.y = unit(0.5,'cm'),
    panel.spacing.x = unit(0.5,'cm'),
    legend.position="right",
    legend.justification="right",
    legend.box.spacing = unit(-0.001, "cm"))
```

# volcano plot for genome-wide RNA-Seq changes TKO/QKO
```r
volc(r=resLFC,x = 'log2FoldChange',y = 'padj',ylab = '-log10(padj)',xlab = "Log2 (FC) gene expression QKO/TKO",ttl = "",cutoff = 2)
```

# correlation heatmaps

## H3K4me1, H3K4me3 and H3K27ac
```bash
multiBigwigSummary bins -b ${samp.bw} -o ${summary}.npz --labels ${sample_label} -bs 500 -bl mm10_blacklist.bed -p 6 -v 
```

## H3K27me3
```bash
multiBigwigSummary bins -b ${samp.bw} -o ${summary}.npz --labels ${sample_label} -bs 10000 -bl mm10_blacklist.bed -p 6 -v 
```

# genome browser track of PA and TKO H3K36me2 using trackplot
```r
bigWigs = c("mMSC_PA_H3K36me2_msNorm.bw",
            "mMSC_TKO_H3K36me2_msNorm.bw")

#Make a table of bigWigs along with ref genome build
bigWigs = read_coldata(bws = bigWigs, build = "mm10")

# peaks
pks=c("mMSC_TKO_strong_enhancers.bed")

#Region to plot
loci = "chr2:37705985-39392970"


t = track_extract(colData = bigWigs, loci = loci,build = "mm10",nthreads = 8,binsize = 50,padding = 10000)

pdf(filename = "mMSC_PA_TKO_enhancer_genome_browser_track.pdf",width = 12,height = 9,res = 600,units = "cm")
track_plot(summary_list = t,col=c("blue","darkkhaki","forestgreen","darkcyan"),track_names = c("PA-H3K36me2","TKO-H3K36me2"),peaks =pks,peaks_track_names = c("Enhancer","Strong enhancer"),draw_gene_track = TRUE,bw_track_height = 2,peaks_track_height = c(1),show_ideogram = FALSE,gene_track_height = 2,left_mar = 3,show_axis = TRUE,gene_fsize = 0.5,y_max = c(0.17,0.17),y_min = c(0,0),collapse_txs = TRUE)
```

# NSD3-KO versus PA comparison
```r
```{r}
load("mMSC_genes_within_large_H3K36me2_domains_in_TKO.RData")
load("mMSC_RNAseq_dds.RData")
# resLFC
dds$condition <- relevel(dds$condition,"PA")
# model.matrix(~ condition
dds <- nbinomWaldTest(dds)
resultsNames <- resultsNames(dds) %>% as.data.frame()
resLFC <- lfcShrink(dds=dds,coef=3,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("gene_id") %>% na.omit()
set.seed(1234)
to_keep=mMSC_genes_within_large_H3K36me2_domains_in_TKO
# genes with K36me2
genes_with_K36me2 <- resLFC %>% .[resLFC$gene_id %in% to_keep$id,] %>% dplyr::select(log2FoldChange) %>% `names<-`(c('with_K36me2_log2FC')) %>% as.data.frame()
nrow(as.data.frame(genes_with_K36me2[genes_with_K36me2$with_K36me2_log2FC > 0,]))
# genes without K36me2
'%!in%' <- function(x,y)!('%in%'(x,y))
genes_without_K36me2 <- resLFC %>% .[resLFC$gene_id %!in% to_keep$id,] %>% dplyr::slice_sample(n=nrow(genes_with_K36me2),replace=FALSE) %>% dplyr::select(log2FoldChange)  %>% `names<-`(c('without_K36me2_log2FC'))
# combine the two
combined_FC <- cbind(genes_with_K36me2,genes_without_K36me2) %>% pivot_longer(cols = c(1:2),names_to=c("type"),values_to=c("log2FoldChange"))

# promoters of genes 
genes <- import.bed("mm10_genes.RData")
promoters <- promoters(genes,upstream=2,downstream=1) %>% as.data.frame()
promoters$gene_id <- gsub("\\..*","",promoters$gene_id)
NSD3_genes <- resLFC %>% .[resLFC$gene_id %in% to_keep$id,]
NSD3_genes_coordinates <- promoters[promoters$gene_id %in% NSD3_genes$gene_id,] %>% mutate(score=0) %>% dplyr::select(c("seqnames","start","end","gene_name","score","strand"))
# full gene bodies
genes <- import.bed("mm10_genes.RData")
gene_bodies <- genes %>% as.data.frame()
gene_bodies$gene_id <- gsub("\\..*","",gene_bodies$gene_id)
NSD3_genes <- resLFC %>% .[resLFC$gene_id %in% to_keep$id,]
NSD3_genes_coordinates <- gene_bodies[gene_bodies$gene_id %in% NSD3_genes$gene_id,] %>% mutate(score=0) %>% dplyr::select(c("seqnames","start","end","gene_name","score","strand"))

size_of_text <- as.numeric(2)
stats <- compare_means(log2FoldChange ~ type,data=combined_FC,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p=",.$p.format))
# run ggplot
ggplot(data = combined_FC,aes(x=type,y=log2FoldChange,fill=type)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1)) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.1) +
  geom_boxplot(width=0.1, color="darkslategray", alpha=1,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2)) +
  scale_fill_manual(values=c("yellow3","blue")) +
  scale_x_discrete(labels=c("with_K36me2_log2FC" = "Genes within\n large H3K36me2\n domains", "without_K36me2_log2FC" = "Genes outside of\n large H3K36me2\n domains")) +
  labs(x="",y="Log2 (FC) gene expression (NSD3KO/PA)") +
  # stats
  geom_signif(textsize = 2,tip_length = 0.004,xmin=stats$group1,xmax=stats$group2,annotations=stats$pvalue,y_position = max(combined_FC$log2FoldChange * 1.05)) +
  annotate("text",y=0.5, x = 1.5, label = paste0(nrow(as.data.frame(genes_with_K36me2[genes_with_K36me2$with_K36me2_log2FC > 0,]))),size=size_of_text) +
  annotate("text",y=0.5, x = 2.5, label = paste0(nrow(as.data.frame(genes_without_K36me2[genes_without_K36me2$without_K36me2_log2FC > 0,]))),size=size_of_text) +
  # count those that go down
   annotate("text",y=-0.5, x = 1.5, label = paste0(nrow(as.data.frame(genes_with_K36me2[genes_with_K36me2$with_K36me2_log2FC < 0,]))),size=size_of_text) +
  annotate("text",y=-0.5, x = 2.5, label = paste0(nrow(as.data.frame(genes_without_K36me2[genes_without_K36me2$without_K36me2_log2FC < 0,]))),size=size_of_text) +
  coord_flip() +
  theme(
    # panel
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y.left  = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.line.x.bottom  = element_line(colour="black"),
    axis.line.x.top  = element_line(colour="black"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
    plot.title = element_text(hjust = 0.5,color = "black",size=6,family="Helvetica"),
    # axis
    axis.title.y= element_text(size=6,family="Helvetica",colour = "black"),
    axis.title.x= element_text(size=6,family="Helvetica",colour = "black"),
    axis.text.x = element_text(size=6,family="Helvetica",colour = "black",hjust = 1,vjust=1),
    axis.text.y=element_text(size=6,family="Helvetica",colour = "black"),
    # labels
    strip.background =element_blank(),
    strip.text.x = element_blank(),
    strip.text = element_text(
      size = 6, color = "black",family = "Helvetica"),
    # legend
    legend.text=element_text(size=6,family = "Helvetica",color = "black"),
    legend.title=element_text(size=6,family="Helvetica",color="black"),
    legend.background = element_rect(fill="white"),
    legend.key=element_rect(fill="white"),
    # margins
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
    panel.spacing = unit(0.5,'cm'),
    panel.spacing.y = unit(0.5,'cm'),
    panel.spacing.x = unit(0.5,'cm'),
    legend.position="right",
    legend.justification="right",
    legend.box.spacing = unit(-0.001, "cm"))
ggsave(filename = "mMSC_log2FC_NSD3KO_PA.pdf",path = "outdir",width = 8,height=5.5,units="cm",dpi = 600,device = "pdf",bg = "white")
```

# agg plot of PA and NSD3-KO on promoters
```bash
computeMatrix reference-point -R promoters_marked_by_H3K36me2_in_TKO.bed \
-S PA_H3K36me2_msNorm.bw NSD3KO_H3K36me2_msNorm.bw \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50 -a 3000 -b 3000 \
-p 16 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
# plot
plotHeight=3.5
plotWidth=3.5
colour_list="blue gold"
plotProfile -m ${mat} -o ${plot}.pdf --dpi 600 --plotHeight ${plotHeight} --plotWidth ${plotWidth} --colors ${colour_list} --refPointLabel "TSS" --yAxisLabel "" --plotFileFormat "pdf" --perGroup --plotTitle "H3K36me2" --regionsLabel " "
```

# agg plot of PA and NSD3-KO on enhancers
```bash
computeMatrix reference-point -R enhancers_marked_by_H3K36me2_in_TKO.bed \
-S PA_H3K36me2_msNorm.bw NSD3KO_H3K36me2_msNorm.bw \
-bl ${blacklist_file} \
--referencePoint center \
--binSize 50 -a 3000 -b 3000 \
-p 16 --samplesLabel ${samp_lab} --verbose \
-o ${mat}
# plot
plotHeight=3.5
plotWidth=3.5
colour_list="blue gold"
plotProfile -m ${mat} -o ${plot}.pdf --dpi 600 --plotHeight ${plotHeight} --plotWidth ${plotWidth} --colors ${colour_list} --refPointLabel "TSS" --yAxisLabel "" --plotFileFormat "pdf" --perGroup --plotTitle "H3K36me2" --regionsLabel " "
```