# H3K9me3 levels within clusB genes

## mMSC
```bash
# run featureCounts to count H3K9me3 signal within genes
featureCounts -a mm10_protein_coding_genes.saf -F SAF -o H3K9me3_mm10_protein_coding_genes.counts -T 6 -p *_H3K9me3.bam
```

```r
rm(list=ls())
'%!in%' <- function(x,y)!('%in%'(x,y))

rawCounts <- fread("H3K9me3_protein_coding_genes.counts")

geneLength <- rawCounts$Length

mat <- rawCounts %>% column_to_rownames("Geneid") %>% dplyr::select(-c(1:5))

metadata <- data.frame(kind=colnames(mat)) %>% mutate(cond = .$kind) %>% column_to_rownames("kind") 
H3K9me3_dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = metadata,
                              design = ~1)

mcols(H3K9me3_dds)$basepairs <- geneLength

mMSC_H3K9me3_FPKM_counts <- DESeq2::fpkm(H3K9me3_dds,robust = FALSE) %>% as.data.frame() %>% tibble::rownames_to_column("symb") %>% `names<-`(c('symb', 'PA',"TKO")) %>% dplyr::mutate(PA = .$PA * WT_H3K9me3_ms) %>% dplyr::mutate(TKO = .$TKO * TKO_H3K9me3_ms) %>% dplyr::mutate(log2FoldChange = log2((.$TKO +1e-3)/(.$PA +1e-3)))

proteinCodingGenes <- fread("mm10_protein_coding_genes.bed")

clustB <- import.bed("mMSC_H3K9me3_clusB.bed")

toFiltGenesWithinClusB <- subsetByOverlaps(proteinCodingGenes,clustB,type = "within") %>% as.data.frame()

clustB_H3K9me3_genes <- mMSC_H3K9me3_FPKM_counts %>% dplyr::filter(ensembl_id %in% toFiltGenesWithinClusB$name) %>% dplyr::mutate(cond="TKO/PA")

fontsize=5
ggplot(data = clustB_H3K9me3_genes,aes(x=cond,y=log2FoldChange,fill=cond)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),linewidth=0.3) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.3,alpha=0.5) +
  geom_boxplot(width=0.15, color="black", alpha=0.1,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("slateblue")) +
  geom_hline(yintercept =0,linetype="dashed",linewidth=0.2) +
  labs(x="",y="H3K9me3 log2 (FC) TKO/PA",title="mMSC genes in\n cluster B regions") +
  geom_vline(xintercept = 2,colour="blue") +
  annotate("text",y=3, x = 1.5, label = paste0(nrow(clustB_H3K9me3_genes[clustB_H3K9me3_genes$log2FoldChange > 0,])),size=1.6) +
   annotate("text",y=-3, x = 1.5, label = paste0(nrow(clustB_H3K9me3_genes[clustB_H3K9me3_genes$log2FoldChange < 0,])),size=1.6) +
  coord_flip() +
  theme_bw()+
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis labels 
      axis.title.y= element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black",hjust = 0.5,vjust=0.5),
      axis.text.y=element_blank(),
      plot.title = element_text(hjust = 0.5,color = "black",size=fontsize,family="Helvetica"),
      # background
      strip.background =element_rect(fill="white"),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
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
      panel.spacing = unit(0.1,'cm'),
      panel.spacing.y = unit(0,'cm'),
      panel.spacing.x = unit(0.1,'cm'))
ggsave(filename = "mMSC_H3K9me3_within_genes_in_clusterB_regions.pdf",path = "outdir",width = 3.1,height=3,units="cm",dpi = 600,device = "pdf")
```

## Cal27
```bash
# run featureCounts to count H3K9me3 signal within genes
featureCounts -a hg38_protein_coding_genes.saf -F SAF -o Cal27_H3K9me3_hg38_protein_coding_genes.counts -T 6 -p *_H3K9me3.bam
```

```r
rm(list=ls())
'%!in%' <- function(x,y)!('%in%'(x,y))

rawCounts <- fread("H3K9me3_protein_coding_genes.counts")

geneLength <- rawCounts$Length

mat <- rawCounts %>% column_to_rownames("Geneid") %>% dplyr::select(-c(1:5))

metadata <- data.frame(kind=colnames(mat)) %>% mutate(cond = .$kind) %>% column_to_rownames("kind") 
H3K9me3_dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = metadata,
                              design = ~1)

mcols(H3K9me3_dds)$basepairs <- geneLength

Cal27_H3K9me3_FPKM_counts <- DESeq2::fpkm(H3K9me3_dds,robust = FALSE) %>% as.data.frame() %>% tibble::rownames_to_column("symb") %>% `names<-`(c('symb', 'PA',"H3K36MOE")) %>% dplyr::mutate(PA = .$PA * WT_H3K9me3_ms) %>% dplyr::mutate(H3K36MOE = .$H3K36MOE * H3K36MOE_H3K9me3_ms) %>% dplyr::mutate(log2FoldChange = log2((.$H3K36MOE +1e-3)/(.$PA +1e-3)))

proteinCodingGenes <- fread("hg38_protein_coding_genes.bed")

clustB <- import.bed("Cal27_H3K9me3_clusB.bed")

toFiltGenesWithinClusB <- subsetByOverlaps(proteinCodingGenes,clustB,type = "within") %>% as.data.frame()

clustB_H3K9me3_genes <- H3K9me3_FPKM_counts %>% dplyr::filter(ensembl_id %in% toFiltGenesWithinClusB$name) %>% dplyr::mutate(cond="H3K36MOE/PA")

fontsize=5
ggplot(data = clustB_H3K9me3_genes,aes(x=cond,y=log2FoldChange,fill=cond)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),linewidth=0.3) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.3,alpha=0.5) +
  geom_boxplot(width=0.15, color="black", alpha=0.1,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("slateblue")) +
  geom_hline(yintercept =0,linetype="dashed",linewidth=0.2) +
  labs(x="",y="H3K9me3 log2 (FC) H3K36MOE/PA",title="Cal27 genes in\n cluster B regions") +
  geom_vline(xintercept = 2,colour="blue") +
  annotate("text",y=3, x = 1.5, label = paste0(nrow(clustB_H3K9me3_genes[clustB_H3K9me3_genes$log2FoldChange > 0,])),size=1.6) +
   annotate("text",y=-3, x = 1.5, label = paste0(nrow(clustB_H3K9me3_genes[clustB_H3K9me3_genes$log2FoldChange < 0,])),size=1.6) +
  coord_flip() +
  theme_bw()+
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis labels 
      axis.title.y= element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black",hjust = 0.5,vjust=0.5),
      axis.text.y=element_blank(),
      plot.title = element_text(hjust = 0.5,color = "black",size=fontsize,family="Helvetica"),
      # background
      strip.background =element_rect(fill="white"),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
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
      panel.spacing = unit(0.1,'cm'),
      panel.spacing.y = unit(0,'cm'),
      panel.spacing.x = unit(0.1,'cm'))
ggsave(filename = "Cal27_H3K9me3_within_genes_in_clusterB_regions.pdf",path = "outdir",width = 3.1,height=3,units="cm",dpi = 600,device = "pdf")
```


# Gene expression changes for genes with significant H3K9me3 loss in their genebodies
```r
'%!in%' <- function(x,y)!('%in%'(x,y))

load("H3K9me3_dds.RData")
# log2FoldChanges
dds_K9me3$condition <- relevel(dds_K9me3$condition,"PA")
dds_K9me3 <- nbinomWaldTest(dds_K9me3)
resultsNames_K9me3 <- resultsNames(dds_K9me3) %>% as.data.frame()
resLFC_K9me3 <- lfcShrink(dds=dds_K9me3,coef=2,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("gene_name") %>% na.omit()
# load RNAseq
load("mMSC_RNASeq_dds.RData")
dds_RNAseq <- dds
rm(dds)
# log2FoldChanges 
dds_RNAseq$condition <- relevel(dds_RNAseq$condition,"NSD_unedit")
dds_RNAseq <- nbinomWaldTest(dds_RNAseq)
resultsNames_RNAseq <- resultsNames(dds_RNAseq) %>% as.data.frame()
# run shrinkage
resLFC_RNAseq <- lfcShrink(dds=dds_RNAseq,coef=4,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("id") %>% na.omit() %>% left_join(.,gene_symbol,by="id")
# filter by FPKM
geneLength <- fread("mm10_geneLengths.tsv")
mcols(dds_RNAseq)$basepairs <- geneLength$length
norm.counts <- DESeq2::fpkm(dds_RNAseq,robust = TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("id") %>% mutate(Mean_exp=(rowMeans(.[grep("TKO|PA", names(.), value = TRUE)]))) %>% dplyr::select(c("id","Mean_exp")) %>% dplyr::filter(Mean_exp > 1)
genes_with_decreased_K9me3 <- resLFC_K9me3 %>% dplyr::filter(log2FoldChange < 0 & padj < 0.05)
# load 
protein_coding_genes_in_clusB <- fread("mm10_protein_coding_genes.bed") %>% unique() %>% as.data.frame() 
colnames(protein_coding_genes_in_clusB)[4] <- "gene_name"
clusB_protein_coding_genes_withLFC <- protein_coding_genes_in_clusB %>% left_join(.,resLFC_RNAseq,by="gene_name") %>% na.omit()

genes_in_clusB <- resLFC_K9me3 %>% .[.$gene_name %in% protein_coding_genes_in_clusB$gene_name,] %>% dplyr::filter(log2FoldChange < 0)
RNAseq_for_genes_that_lost_K9me3 <- resLFC_RNAseq %>% .[.$gene_name %in% genes_with_decreased_K9me3$gene_name,]%>% .[.$gene_name %in% protein_coding_genes_in_clusB$gene_name,] %>% dplyr::select(log2FoldChange)
RNAseq_for_genes_that_lost_K9me3$cond <- "Genes with significant \nloss of H3K9me3\nin their gene body"
# random genes
random_genes <- resLFC_K9me3 %>% .[.$gene_name %!in% genes_with_decreased_K9me3$gene_name,]
# random genes
set.seed(1)
RNAseq_for_random_genes <- resLFC_RNAseq %>% .[.$gene_name %in% random_genes$gene_name,] %>% dplyr::select(log2FoldChange) %>% dplyr::slice_sample(n=nrow(RNAseq_for_genes_that_lost_K9me3),replace=FALSE)
RNAseq_for_random_genes$cond <- "Random genes"
# combine
combined_FC <- rbind(RNAseq_for_genes_that_lost_K9me3,RNAseq_for_random_genes)
combined_FC$cond <- factor(combined_FC$cond,levels = c("Random genes","Genes with significant \nloss of H3K9me3\nin their gene body"))

yaxis_label <- "Gene expression \nLog2 (FC) TKO/PA"
fontsize=5
labelsize=1.8
# counting the number that go up and down
second_group <- RNAseq_for_genes_that_lost_K9me3
first_group <- RNAseq_for_random_genes
# run stats
stats <- compare_means(log2FoldChange ~ cond,data=combined_FC,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p=",.$p.format))
# run ggplot
ggplot(data = combined_FC,aes(x=cond,y=log2FoldChange,fill=cond)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),linewidth=0.1) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.1,alpha=0.5) +
  geom_boxplot(width=0.15, color="black", alpha=0.1,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("lightcoral","slateblue")) +
  geom_hline(yintercept =0,linetype="dashed",linewidth=0.1) +
  geom_hline(yintercept =16,color="white") +
  labs(x="",y=yaxis_label) +
  # stats
  geom_signif(textsize = 1.8,tip_length = 0.03,xmin=stats$group1,xmax=stats$group2,annotations=stats$p.signif,y_position = max(c(combined_FC$log2FoldChange * 1.15,combined_FC$log2FoldChange * 1.15,combined_FC$log2FoldChange*1.15)),size=0.1) +
  annotate("text",y=2, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange > 0,])),size=labelsize) +
  annotate("text",y=2, x = 2.5, label = paste0(nrow(second_group[second_group$log2FoldChange > 0,])),size=labelsize) +
  annotate("text",y=-2, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange < 0,])),size=labelsize) +
  annotate("text",y=-2, x = 2.5, label = paste0(nrow(second_group[second_group$log2FoldChange < 0,])),size=labelsize) +
  coord_flip() +
    theme(
      # panel
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis labels 
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black",hjust = 0.5,vjust=0.5),
      axis.text.y=element_text(size=fontsize-1,family="Helvetica",colour = "black"),
      plot.title = element_text(hjust = 0.5,color = "black",size=fontsize,family="Helvetica"),
      axis.ticks = element_line(size=0.1),
      # background
      strip.background =element_rect(fill="white"),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
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
      panel.spacing = unit(0.1,'cm'),
      panel.spacing.y = unit(0,'cm'),
      panel.spacing.x = unit(0.1,'cm'))

ggsave(filename = "TKO_PA.significant_loss_of_K9me3.log2FoldChange_gene_expression.pdf",path="outdir",device = "pdf",units = "cm",width = 5,height=3.8,dpi = 600,bg="white")
```

# overlap for genes losing H3K9me3
```r
mMSC_protein_coding_genes_in_clusB <- fread("exclusive_protein_coding_genes_within_clusB.bed") %>% unique() %>% `names<-`(c('chr', 'start','end','gene_symbol','score','strand'))

# generate mouse-human protein coding orthologs
mm10_protein_coding_genes <- fread("mm10.protein_coding_genes.with_gene_symbols.csv")

mouse_human_orthologs <- fread("mouse_human_orthologs.txt") %>% dplyr::select(c(`Gene stable ID`,`Mouse gene stable ID`)) %>% `names<-`(c('human_ensembl_id', 'id'))

mm10_protein_coding_genes_with_human_orthologs <- mm10_protein_coding_genes %>% left_join(mouse_human_orthologs,by="id") %>% na.omit()

load("mMSC_genesLosingK9me3_inClusB.RData")
mMSC_genes_losingK9me3 <- K9me3_for_genes_in_clusB %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::select(gene_name) %>% mutate(gene_symbol = .$gene_name) %>% dplyr::filter(gene_name %in% mm10_protein_coding_genes$gene_symbol) %>% dplyr::select(c("gene_symbol")) %>% left_join(.,mm10_protein_coding_genes,keep = FALSE,unmatched = "drop",relationship = "one-to-one",multiple = "first") %>% left_join(.,mouse_human_orthologs,keep = FALSE,unmatched = "drop",relationship = "one-to-one",multiple = "first")

load("Cal27_clusB_H3K9me3_genes.RData")
hg38_gene_symbols <- fread("~/Documents/HNSCC_K36me2/ref/hg38_protein_coding_genes.bed") %>% dplyr::select(c("gene_id","gene_name"))
Cal27_genes_losingK9me3 <- clustB_H3K9me3_genes %>% dplyr::filter(log2FoldChange < 0)

overlap <- intersect(Cal27_genes_losingK9me3$ensembl_id,mMSC_genes_losingK9me3$human_ensembl_id) %>% as.data.frame() %>% `names<-`(c('gene_id')) %>% left_join(hg38_gene_symbols,by="gene_id")
```


# TE transcriptional changes in clusB regions
```bash
# run featureCounts on individual TEs
featureCounts -a individual_TE.saf -F SAF -s 2 -o individual_TE.counts -T 6 -p PA_RNAseq_rep1.bam PA_RNAseq_rep2.bam PA_RNAseq_rep3.bam TKO_RNAseq_rep1.bam TKO_RNAseq_rep2.bam TKO_RNAseq_rep3.bam
```

```r
counts <- fread("individual_TE.counts")
df <- counts %>% dplyr::select(1,7:12) %>% column_to_rownames("Geneid")
mat <- df %>% as.matrix()
metadata <- data.frame(kind=colnames(mat))
metadata$condition <- metadata$kind
metadata$condition <- factor(metadata$condition,levels=c("PA","TKO"))
# run DESeq2
metadata <- tibble::column_to_rownames(metadata,"kind")
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = metadata,
                              design = ~condition)
dds <- DESeq(dds)


rmsk <- fread("mm10_UCSC_repeatMasker.txt") %>%
  subset(.,(repClass %in% c("SINE","LINE","LTR","DNA","Retroposon")))
# chrom sizes
keep <- fread("mm10.chrom.sizes") %>%
  setNames(c("chr","seqlength"))
gn <- keep %>% {Seqinfo(.$chr, .$seqlength)}
rmsk_sub <- rmsk[,c("repClass","repFamily","repName","genoName","genoStart","genoEnd","strand")] %>% as.data.frame() %>% subset(.,(genoName %in% keep$chr))
rmsk_sub$TE_id <- c(paste0("TE_",1:(nrow(rmsk_sub))))


resLFC_individual_TE <- lfcShrink(dds=dds,coef=2,type ="apeglm") %>% as.data.frame() %>% rownames_to_column("gene_id") %>% na.omit() %>% left_join(rmsk_sub,by="TE_id") %>% makeGRangesFromDataFrame(df = .,keep.extra.columns = TRUE,seqnames.field = "genoName",start.field = "genoStart",end.field = "genoEnd",ignore.strand = TRUE)

clusA <- import.bed("mMSC_H3K9me3_clusA.bed")
clusB <- import.bed("mMSC_H3K9me3_clusB.bed")

TEs_within_clusB <- resLFC_individual_TE[overlapsAny(resLFC_individual_TE,clusB,type = "within") & !overlapsAny(resLFC_individual_TE, clusA),] %>% as.data.frame()

volc <- function(r, x, y, ylab,xlab, ttl,cutoff,fontsize=5) {
  d <- as.data.frame(r) %>%
    dplyr::rename(x = !!x, y = !!y) %>% dplyr::filter(!is.na(x) & !is.na(y)) %>%
    mutate(kind = case_when((x > cutoff & y < .05) ~ paste('FDR<0.05\nlog2FC>',cutoff,sep = "",collapse = ""),
                            (x < -cutoff & y < .05) ~ paste('FDR<0.05\nlog2FC<',-cutoff,sep = "",collapse = ""),(cutoff > abs(x) & y < .05) ~ ('FDR<0.05'),
                            y > .05 ~ 'NS'),
           y = -log10(y))
  d$kind <- factor(d$kind, levels = c(paste('FDR<0.05\nlog2FC>',cutoff,sep = "",collapse = ""),paste('FDR<0.05\nlog2FC<',-cutoff,sep = "",collapse = ""),'FDR<0.05','NS'))
  ct_up <- d %>% 
    dplyr::filter(kind == paste('FDR<0.05\nlog2FC>',cutoff,sep = "",collapse = "")) %>%
    mutate(up = x > 0) %>%
    dplyr::count(up) %>%
    mutate(x = ifelse(up, Inf, -Inf),
           y = Inf,
           h = as.numeric(up))
  ct_down <- d %>% 
    dplyr::filter(kind == paste('FDR<0.05\nlog2FC<',-cutoff,sep = "",collapse = "")) %>%
    mutate(up = x > 0) %>%
    dplyr::count(up) %>%
    mutate(x = ifelse(up, Inf, -Inf),
           y = Inf,
           h = as.numeric(up))
  ggplot(d, aes(x, y, color = kind)) +
    geom_vline(xintercept = c(-cutoff, cutoff), linetype = 'dashed',linewidth=0.2) +
    geom_hline(yintercept = -log10(0.05),linetype = 'dashed',linewidth=0.2) +
    geom_point(size=0.1,alpha=0.6) +
    geom_label(aes(x = x, y = y, label = n, hjust = h),
              vjust = 1, data = ct_up, inherit.aes = F,size=1.7) + scale_y_continuous() + geom_label(aes(x = x, y = y, label = n, hjust = h),
              vjust = 1, data = ct_down, inherit.aes = F,size=1.7) +
    scale_color_manual(values = c('purple','darkorange4','forestgreen','darkgrey')) +
    labs(x = xlab, y = ylab, title = ttl) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=fontsize,family="Helvetica",colour = "black"),
    axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
    axis.text.x=element_text(size=fontsize,family="Helvetica",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "black",size=fontsize+1,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = fontsize-1, color = "black",family = "Helvetica"),
    # legend
    legend.text=element_text(size=4.5,family = "Helvetica",color = "black"),
    legend.title=element_blank(),
    legend.background = element_rect(fill="white",size = 4),
    legend.key=element_rect(fill="white",size = 2),
    legend.margin = unit(c(0,0,0,0),"mm"),
    legend.spacing = unit(c(0,0,0,0),"mm"),
    legend.position = "top",
    legend.justification = "right",
    legend.box.margin = margin(-5,-5,-5,-5),
    legend.key.size = unit(0,"mm"),
    plot.margin = unit(c(2.4,2.4,2.4,2.4), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0.1,'cm'),
    panel.spacing.x = unit(0.1,'cm'))
}
volc(r=TEs_within_clusB,x = 'log2FoldChange',y = 'padj',ylab = '-log10(padj)',xlab = "Log2 (FC) Expression TKO/PA",ttl = "TEs in cluster B regions",cutoff = 1)
ggsave(filename = "TKO_PA.individual_TEs_in_clusB.volcano_plot.pdf",path="outdir",device = "pdf",units = "cm",width = 4.5,height=4.7,dpi = 600,bg="white")
```

# subfamily TE analysis
```r
load("annotated_rmsk.RData")
colnames(complete_annotated_rmsk)[8] <- "gene_id"
load("resLFC_subfamily.RData")
# resLFC subfamilies genome-wide
res_subfamily_TE <- resLFC_subfamilies %>% dplyr::select(c("TE_subfamily","log2FoldChange","padj")) %>% `names<-`(c('repName', 'subFamily_log2FC','subFamily_padj'))
# load clusB
clusB <- import.bed("mMSC_clusB.bed")
# left bind individual TEs within clusB with subfamily
TEs_within_clusB_with_subfamily_info <- TEs_within_clusB %>% left_join(complete_annotated_rmsk,by="gene_id") %>% left_join(res_subfamily_TE,by="repName")
clusB_subfamilies <- TEs_within_clusB_with_subfamily_info %>% dplyr::select(c("gene_id","repFamily","repName","repClass","subFamily_log2FC","subFamily_padj")) %>% .[!duplicated(.[,c("repName")]),]
colnames(clusB_subfamilies)[3] <- "TE_subfamily"

volc <- function(r, x, y, ylab,xlab, ttl,cutoff) {
  d <- as.data.frame(r) %>%
    dplyr::rename(x = !!x, y = !!y) %>% dplyr::filter(!is.na(x) & !is.na(y)) %>%
    mutate(kind = case_when((x > cutoff & y < .05) ~ paste('FDR<0.05,log2FC>',cutoff,sep = "",collapse = ""),
                            (x < -cutoff & y < .05) ~ paste('FDR<0.05,log2FC<',-cutoff,sep = "",collapse = ""),(cutoff > abs(x) & y < .05) ~ ('FDR<0.05'),
                            y > .05 ~ 'NS'),
           y = -log10(y))
  d$kind <- factor(d$kind, levels = c(paste('FDR<0.05,log2FC>',cutoff,sep = "",collapse = ""),paste('FDR<0.05,log2FC<',-cutoff,sep = "",collapse = ""),'FDR<0.05','NS'))
  ct_up <- d %>% 
    dplyr::filter(kind == paste('FDR<0.05,log2FC>',cutoff,sep = "",collapse = "")) %>%
    mutate(up = x > 0) %>%
    dplyr::count(up) %>%
    mutate(x = ifelse(up, Inf, -Inf),
           y = Inf,
           h = as.numeric(up))
  ct_down <- d %>% 
    dplyr::filter(kind == paste('FDR<0.05,log2FC<',-cutoff,sep = "",collapse = "")) %>%
    mutate(up = x > 0) %>%
    dplyr::count(up) %>%
    mutate(x = ifelse(up, Inf, -Inf),
           y = Inf,
           h = as.numeric(up))
  ggplot(d, aes(x, y, color = kind)) +
    geom_vline(xintercept = c(-cutoff, cutoff), linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05),linetype = 'dashed') +
    geom_point(size=0.2) +
    geom_label(aes(x = x, y = y, label = n, hjust = h),
              vjust = 1, data = ct_up, inherit.aes = F,size=2.5) + scale_y_continuous() + geom_label(aes(x = x, y = y, label = n, hjust = h),
              vjust = 1, data = ct_down, inherit.aes = F,size=2.5) +
    geom_text_repel(aes(label = TE_subfamily), data = d %>% dplyr::filter((kind %in% c(paste('FDR<0.05,log2FC>',cutoff,sep = "",collapse = ""),paste('FDR<0.05,log2FC<',-cutoff,sep = "",collapse = "")))),max.overlaps = 5,show.legend = F, min.segment.length = 0,size=1.5) +
    scale_color_manual(values = c('purple','darkorange4','forestgreen','darkgrey')) +
    labs(x = xlab, y = ylab, title = ttl) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=6,family="Helvetica",colour = "black"),
    axis.text.y=element_text(size=6,family="Helvetica",colour = "black"),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "black",size=6,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = 8, color = "black",family = "Helvetica"),
    # legend
    legend.text=element_text(size=4.5,family = "Helvetica",color = "black"),
    legend.title=element_blank(),
    legend.background = element_rect(fill="white",size = 4),
    legend.key=element_rect(fill="white",size = 4),
    legend.margin = unit(c(0,0,0,0),"mm"),
    legend.spacing = unit(c(0,0,0,0),"mm"),
    legend.position = "top",
    legend.justification = "right",
    legend.box.margin = margin(-5,-5,-5,-5),
    legend.key.size = unit(0,"mm"),
    plot.margin = unit(c(2.4,2.4,2.4,2.4), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0.1,'cm'),
    panel.spacing.x = unit(0.1,'cm'))
}

volc(r=clusB_subfamilies,x = 'subFamily_log2FC',y = 'subFamily_padj',ylab = '-log10(padj)',xlab = "Log2 (FC) expression TKO/PA",ttl = "TE subfamilies in cluster B regions",cutoff = 1)
ggsave(filename = "TKO_PA.subfamily_TEs_in_clusB.volcano_plot.pdf",path="outdir",device = "pdf",units = "cm",width = 5.5,height=5,dpi = 600,bg="white")

res <- results(dds,contrast=c("condition","TKO","PA")) %>% as.data.frame() %>% na.omit() %>% rownames_to_column("repName")
TE_family <- rmsk %>% dplyr::select(c("repName","repFamily"))
TE_family <- TE_family[TE_family$repName %in% res$repName,] %>% unique()
res <- res %>% left_join(.,TE_family,by="repName")
upregs <- res %>% dplyr::filter(res$log2FoldChange > 1 & res$padj < 0.05) %>% dplyr::group_by(repFamily) %>% dplyr::count(repFamily) %>% dplyr::mutate(Expression = "pos") %>% dplyr::filter(n>=3)
downregs <- res %>% dplyr::filter(res$log2FoldChange < -1 & res$padj < 0.05) %>% dplyr::group_by(repFamily) %>% dplyr::count(repFamily) %>% dplyr::mutate(Expression = "neg") %>% dplyr::filter(n>=3)
combined_regs <- bind_rows(upregs,downregs) %>% dplyr::group_by(Expression)
combined_regs$Expression <- factor(combined_regs$Expression,levels=c("neg","pos"))
# cluster B
upregs_clusB <- clusB_subfamilies %>% dplyr::filter(clusB_subfamilies$subFamily_log2FC > 1 & clusB_subfamilies$subFamily_padj < 0.05) %>% dplyr::group_by(repFamily) %>% dplyr::count(repFamily) %>% dplyr::mutate(Expression = "pos") %>% dplyr::filter(n>=3)
downregs_clusB <- clusB_subfamilies  %>% dplyr::filter(clusB_subfamilies$subFamily_log2FC  < -1 & clusB_subfamilies$subFamily_padj< 0.05) %>% dplyr::group_by(repFamily) %>% dplyr::count(repFamily) %>% dplyr::mutate(Expression = "neg") %>% dplyr::filter(n>=3)
combined_regs_clusB <- bind_rows(upregs_clusB,downregs_clusB) %>% dplyr::group_by(Expression)
combined_regs_clusB$Expression <- factor(combined_regs_clusB$Expression,levels=c("neg","pos"))

ggplot(data=combined_regs_clusB,aes(x=factor(repFamily),y=n,fill=Expression)) +
  geom_col(position = position_dodge2(width=0.9,preserve="single"))+
  scale_fill_manual(values=c("purple"),labels=c("Upregulated"),name="Expression") +
  labs(x="Family",y="Count",title="Subfamilies significantly regulated\nin one direction in cluster B") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=6,family="Helvetica",colour = "black"),
    axis.text.y=element_text(size=6,family="Helvetica",colour = "black"),
    axis.text.x=element_text(size=5,family="Helvetica",colour = "black",angle=90,hjust=1,vjust=0.3),
    axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
    plot.title = element_text(hjust = 0.5,color = "black",size=6,family="Helvetica"),
    strip.background =element_rect(fill="white"),
    strip.text = element_text(
        size = 8, color = "black",family = "Helvetica"),
    # legend
    legend.text=element_text(size=4.5,family = "Helvetica",color = "black"),
    legend.title=element_text(size=4.5,family = "Helvetica",color = "black"),
    legend.margin = unit(c(0,0,0,0),"mm"),
    legend.spacing = unit(c(0,0,0,0),"mm"),
    legend.key.size = unit(2,"mm"),
    legend.position = "top",
    legend.justification = "right",
    legend.box.margin = margin(-5,-5,-5,-5),
    plot.margin = unit(c(2.4,2.4,2.4,2.4), "mm"),
    panel.spacing = unit(0.1,'cm'),
    panel.spacing.y = unit(0.1,'cm'),
    panel.spacing.x = unit(0.1,'cm'))
ggsave(filename = "mMSC_TKO_PA_subfamilies_regulated_onedirection_clusterB.pdf",path="outdir",device = "pdf",units = "cm",width = 5,height=6,dpi = 600,bg="white")
```


# run over-representation pathway analysis using ClusterProfiler
```r
hg38_protein_coding_genes <- fread("hg38_protein_coding_genes.bed")

overlap <- intersect(Cal27_genes_in_clusB,mMSC_orthologousGenes_losingK9me3) %>% as.data.frame() 

query_ids <- mapIds(org.Hs.eg.db, keys = overlap$gene_id, keytype="ENSEMBL", column = "ENTREZID") %>% as.data.frame()
bg_ids <- mapIds(org.Hs.eg.db, keys = hg38_gene_symbols$gene_id, keytype="ENSEMBL", column = "ENTREZID") %>% as.data.frame()

ekegg <- enrichKEGG(gene=query_ids,organism='hsa',pvalueCutoff=0.05,universe = bg_ids)

ego <- enrichGO(gene = query_ids, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "ALL", 
                    pAdjustMethod = "fdr", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE,universe = bg_ids)
```

# run partial correlation network analysis
```r
# in mMSC
resLFC_RNAseq <- resLFC_RNAseq %>% dplyr::filter(padj < 0.05)
resLFC_H3K27me3 <- resLFC_H3K27me3 %>% dplyr::filter(id %in% resLFC_RNAseq$id & id %in% resLFC_H3K9me3$id) %>% mutate(H3K27me3 = log2FoldChange)
resLFC_H3K9me3 <- resLFC_H3K9me3 %>% dplyr::filter(id %in% resLFC_RNAseq$id & id %in% resLFC_H3K27me3$id) %>% mutate(H3K9me3 = log2FoldChange)
resLFC_RNAseq <- resLFC_RNAseq %>% dplyr::filter(id %in% resLFC_H3K27me3$id & id %in% resLFC_H3K9me3$id) %>% mutate(RNAseq_log2FC = log2FoldChange)

# agg
agg <- left_join(resLFC_H3K27me3 %>% dplyr::select(c("id","H3K27me3")),resLFC_RNAseq %>% dplyr::select(c("id","RNAseq_log2FC"))) %>% left_join(resLFC_H3K9me3 %>% dplyr::select(c("id","H3K9me3")))

agg$H3K27me3_centered <- scale(agg$H3K27me3,center=TRUE,scale=FALSE)
agg$H3K9me3_centered <- scale(agg$H3K9me3,center=TRUE,scale=FALSE)

agg_mat <- agg %>% dplyr::select(c("RNAseq_log2FC","H3K27me3","H3K9me3")) %>% as.data.frame()

cor <- agg_mat %>%
  set_names("FC","K27","K9") %>%
  cor_auto() %>% 
  qgraph(graph = "cor", threshold = "fdr",sampleSize = nrow(agg_mat), alpha = 0.05, layout = "spring", title = "Partial correlation network", details = TRUE, edge.labels = TRUE)
```

# Cal27 gene expression changes
```r
hg38ProteinCodingGenes <- import.bed("hg38_protein_coding_genes.bed") %>% as.data.frame() %>% mutate(start = .$start - 1)

raw_counts <- fread("Cal27_RNAseq.counts") %>% dplyr::filter(Geneid %in% hg38ProteinCodingGenes$name)

gene_length <- raw_counts$Length

mat <- raw_counts %>%
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix() 

metadata <- data.frame(kind = colnames(mat))

metadata$condition <- gsub(pattern = "_[0-9]", "", metadata$kind)
# factorization
metadata$condition <- factor(metadata$condition, levels = c("WT","H3K36M_OE"))
# run DESeq2
metadata <- tibble::column_to_rownames(metadata, "kind")
dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData = metadata,
  design = ~condition
)
dds <- DESeq(dds)

resNames <- resultsNames(object = dds) %>% as.data.frame()

Cal27_PA_K36M_resLFC <- lfcShrink(dds = dds,coef = 2,type = "apeglm") %>% as.data.frame() %>% na.omit() %>% rownames_to_column("Geneid")

hg38_symb <- fread("~/Documents/HNSCC_K36me2/ref/hg19_gene_symbol_length.tsv") %>% dplyr::select(c(1:2)) %>% `names<-`(c('symb', 'Geneid'))

colnames(hg38ProteinCodingGenes)[6] <- "Geneid"

Cal27_PA_K36M_geneExpressionChanges <- hg38ProteinCodingGenes %>% dplyr::left_join(.,Cal27_PA_K36M_resLFC %>% dplyr::select(c("Geneid","log2FoldChange"))) %>% na.omit() %>% dplyr::select(c("seqnames","start","end","Geneid","log2FoldChange","strand")) %>% dplyr::left_join(hg38_symb,by="Geneid")
```

```bash
bedtools intersect -a Cal27_PA_K36M_geneExpressionChanges -b Cal27_clusB.bed -wa -f 1.0 > Cal27_genes_within_clusB.bed
```

```r
set.seed(46)
Cal27_genesWithinClusB <- fread("Cal27_genes_within_clusB.bed") %>% `names<-`(c("chr","start","end","ensembl_id","log2FoldChange","strand","symb")) 

RandomGenes <- Cal27_PA_K36M_resLFC %>% dplyr::filter(!Geneid %in% Cal27_genesWithinClusB$ensembl_id) %>% slice_sample(n=nrow(Cal27_genesWithinClusB),replace=FALSE)

agg <- rbind(Cal27_genesWithinClusB %>% dplyr::select(log2FoldChange) %>% mutate(cond="GenesWithinClusB"),RandomGenes %>% dplyr::select(log2FoldChange)%>%mutate(cond="RandomGenes"))

agg$cond <- factor(agg$cond,levels = c("RandomGenes","GenesWithinClusB"))

fontsize=4.5
second_group <- agg %>% dplyr::filter(cond == "GenesWithinClusB")
first_group <- agg %>% dplyr::filter(cond == "RandomGenes")
# run stats
stats <- compare_means(log2FoldChange ~ cond,data=agg,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p=",.$p.format))
# run ggplot
ggplot(data = agg,aes(x=cond,y=log2FoldChange,fill=cond)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),linewidth=0.3) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.1,alpha=0.7) +
  geom_boxplot(width=0.2, color="black", alpha=0.1,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("slateblue","lightcoral")) +
  geom_hline(yintercept =0,linetype="dashed",linewidth=0.1) +
  labs(x="",y="Gene expression\nlog2 (FC) H3K36M-OE/PA",title="Cal27") +
  # stats
  geom_signif(textsize = 2,tip_length = 0.01,xmin=stats$group1,xmax=stats$group2,annotations=stats$p.signif,y_position = max(c(agg$log2FoldChange * 1.05,agg$log2FoldChange * 1.05,agg$log2FoldChange*1.05)),size=0.1) +
  # count the number that go up
  annotate("text",y=2, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange > 0,])),size=2) +
  annotate("text",y=2, x = 2.5, label = paste0(nrow(second_group[second_group$log2FoldChange > 0,])),size=2) +
  # count those that go down
   annotate("text",y=-2, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange < 0,])),size=2) +
  annotate("text",y=-2, x = 2.5, label = paste0(nrow(second_group[second_group$log2FoldChange < 0,])),size=2) +
  coord_flip() +
  geom_hline(yintercept = 11,colour="white")+
  scale_x_discrete(name = "",
                   labels=c("GenesWithinClusB" = "Genes within\nclusB",
                            "RandomGenes" = "Random\ngenes"))+
    theme(
      # panel
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis labels 
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black",hjust = 0.5,vjust=0.5),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
      plot.title = element_text(hjust = 0.5,color = "black",size=7,family="Helvetica"),
      # background
      strip.background =element_rect(fill="white"),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
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
      panel.spacing = unit(0.1,'cm'),
      panel.spacing.y = unit(0,'cm'),
      panel.spacing.x = unit(0.1,'cm'))
ggsave(filename = "Cal27_genes_losing_K9me3_log2FC_gene_expression_comparison.pdf",path="outdir",device = "pdf",units = "cm",width = 4,height=3.5,dpi = 600,bg="white")
```

# Detroit gene expression changes
```r

hg38ProteinCodingGenes <- import.bed("hg38_protein_coding_genes.bed") %>% as.data.frame() %>% mutate(start = .$start - 1)

raw_counts <- fread("Detroit_RNAseq.counts") %>% dplyr::filter(Geneid %in% hg38ProteinCodingGenes$name)

gene_length <- raw_counts$Length

mat <- raw_counts %>%
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix() 

metadata <- data.frame(kind = colnames(mat))

metadata$condition <- gsub(pattern = "_[0-9]", "", metadata$kind)
# factorization
metadata$condition <- factor(metadata$condition, levels = c("WT","H3K36M_OE"))
# run DESeq2
metadata <- tibble::column_to_rownames(metadata, "kind")
dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData = metadata,
  design = ~condition
)
dds <- DESeq(dds)

resNames <- resultsNames(object = dds) %>% as.data.frame()

Detroit_PA_K36M_resLFC <- lfcShrink(dds = dds,coef = 2,type = "apeglm") %>% as.data.frame() %>% na.omit() %>% rownames_to_column("Geneid")

hg38_symb <- fread("~/Documents/HNSCC_K36me2/ref/hg19_gene_symbol_length.tsv") %>% dplyr::select(c(1:2)) %>% `names<-`(c('symb', 'Geneid'))

colnames(hg38ProteinCodingGenes)[6] <- "Geneid"

Detroit_PA_K36M_geneExpressionChanges <- hg38ProteinCodingGenes %>% dplyr::left_join(.,Detroit_PA_K36M_resLFC %>% dplyr::select(c("Geneid","log2FoldChange"))) %>% na.omit() %>% dplyr::select(c("seqnames","start","end","Geneid","log2FoldChange","strand")) %>% dplyr::left_join(hg38_symb,by="Geneid")
```

```bash
bedtools intersect -a Detroit_PA_K36M_geneExpressionChanges -b Detroit_clusB.bed -wa -f 1.0 > Detroit_genes_within_clusB.bed
```

```r
set.seed(46)
Detroit_genesWithinClusB <- fread("Detroit_genes_within_clusB.bed") %>% `names<-`(c("chr","start","end","ensembl_id","log2FoldChange","strand","symb")) 

RandomGenes <- Detroit_PA_K36M_resLFC %>% dplyr::filter(!Geneid %in% Detroit_genesWithinClusB$ensembl_id) %>% slice_sample(n=nrow(Detroit_genesWithinClusB),replace=FALSE)

agg <- rbind(Detroit_genesWithinClusB %>% dplyr::select(log2FoldChange) %>% mutate(cond="GenesWithinClusB"),RandomGenes %>% dplyr::select(log2FoldChange)%>%mutate(cond="RandomGenes"))

agg$cond <- factor(agg$cond,levels = c("RandomGenes","GenesWithinClusB"))

fontsize=4.5
second_group <- agg %>% dplyr::filter(cond == "GenesWithinClusB")
first_group <- agg %>% dplyr::filter(cond == "RandomGenes")
# run stats
stats <- compare_means(log2FoldChange ~ cond,data=agg,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p=",.$p.format))
# run ggplot
ggplot(data = agg,aes(x=cond,y=log2FoldChange,fill=cond)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),linewidth=0.3) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.1,alpha=0.7) +
  geom_boxplot(width=0.2, color="black", alpha=0.1,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("slateblue","lightcoral")) +
  geom_hline(yintercept =0,linetype="dashed",linewidth=0.1) +
  labs(x="",y="Gene expression\nlog2 (FC) H3K36M-OE/PA",title="Detroit") +
  # stats
  geom_signif(textsize = 2,tip_length = 0.01,xmin=stats$group1,xmax=stats$group2,annotations=stats$p.signif,y_position = max(c(agg$log2FoldChange * 1.05,agg$log2FoldChange * 1.05,agg$log2FoldChange*1.05)),size=0.1) +
  # count the number that go up
  annotate("text",y=2, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange > 0,])),size=2) +
  annotate("text",y=2, x = 2.5, label = paste0(nrow(second_group[second_group$log2FoldChange > 0,])),size=2) +
  # count those that go down
   annotate("text",y=-2, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange < 0,])),size=2) +
  annotate("text",y=-2, x = 2.5, label = paste0(nrow(second_group[second_group$log2FoldChange < 0,])),size=2) +
  coord_flip() +
  geom_hline(yintercept = 11,colour="white")+
  scale_x_discrete(name = "",
                   labels=c("GenesWithinClusB" = "Genes within\nclusB",
                            "RandomGenes" = "Random\ngenes"))+
    theme(
      # panel
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # axis labels 
      axis.title.y= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.title.x= element_text(size=fontsize,family="Helvetica",colour = "black"),
      axis.text.x = element_text(size=fontsize,family="Helvetica",colour = "black",hjust = 0.5,vjust=0.5),
      axis.text.y=element_text(size=fontsize,family="Helvetica",colour = "black"),
      plot.title = element_text(hjust = 0.5,color = "black",size=7,family="Helvetica"),
      # background
      strip.background =element_rect(fill="white"),
      strip.text = element_text(
        size = fontsize, color = "black",family = "Helvetica"),
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
      panel.spacing = unit(0.1,'cm'),
      panel.spacing.y = unit(0,'cm'),
      panel.spacing.x = unit(0.1,'cm'))
ggsave(filename = "Detroit_genes_losing_K9me3_log2FC_gene_expression_comparison.pdf",path="outdir",device = "pdf",units = "cm",width = 4,height=3.5,dpi = 600,bg="white")
```