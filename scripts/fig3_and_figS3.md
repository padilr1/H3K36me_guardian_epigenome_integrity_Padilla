# split genes into quartiles
```r
# normalize using FPKM
load("mMSC_RNASeq_dds.RData")
mcols(dds)$basepairs <- length$length
FPKMcounts <- DESeq2::fpkm(dds,robust = TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("id") %>% mutate(PA_mean=(rowMeans(.[grep("PA", names(.), value = TRUE)])))

load("mm10_genes.RData")
genes.df <- genes %>% as.data.frame()
genes.df$gene_id <- gsub("[.][0-9]*","",genes.df$gene_id)
gene_coordinates <- genes.df %>% dplyr::filter(gene_type=="protein_coding") %>% dplyr::select(c("seqnames","start","end","width","strand","gene_id"))

df <- finalset_genes %>% as.data.frame()
# Convert Chr and Strand columns to characters
df$chr <- as.character(df$chr)
df$strand <- as.character(df$strnd)
df$start <- as.numeric(df$start)
df$end <- as.numeric(df$end)

# Function to remove adjacent genes within a 2kb range on the same chr and strand
remove_adjacent_genes <- function(df) {
  # Sort dataframe by chromosome, strand, and start coordinates
  sorted_indices <- order(df$chr, df$strand, df$start)
  df <- df[sorted_indices, ]
  
  # Initialize an empty vector to store indices of genes to be removed
  to_remove <- c()
  
  # Iterate through the dataframe to identify adjacent genes
  for (i in 1:(nrow(df) - 1)) {
    # Check if the genes are on the same chromosome and strand
    if (df$chr[i] == df$chr[i + 1] && df$strand[i] == df$strand[i + 1]) {
      # Check if the distance between the end of the current gene and the start of the next gene is within 2kb
      if ((df$start[i + 1] - df$end[i]) <= 2000) {
        # Mark both genes for removal
        to_remove <- c(to_remove, i, i + 1)
      }
    }
  }
  
  # Remove duplicate indices and remove adjacent genes from the dataframe
  df <- df[-unique(to_remove), ]
  
  return(df)
}

# Call the function to remove adjacent genes
result_df <- remove_adjacent_genes(df)

# median size of genes is 18785
# 25th is 6854
# 75th is 48614
quantile(result_df$width)
final_df <- result_df %>% dplyr::filter(width >= 6854 & width <= 48614)

# stratifying the genes
zero <-  final_df %>%
  dplyr::filter(PA_mean == 0)
exp_stratified <- final_df %>%
  dplyr::filter(PA_mean > 0.1) %>%
  mutate(quantile = ntile(PA_mean,3))

# split genes
zero_expression <- zero 
lowly_expressed <- exp_stratified %>% dplyr::filter(quantile == 1)
highly_expressed <- exp_stratified %>% dplyr::filter(quantile != 1) %>% dplyr::arrange(desc(PA_mean)) %>% dplyr::slice(1:2417)
# reorder columns
silent_genes <- zero_expression[,c("chr","start","end","gene_id","PA_mean","strand","symbol")] %>% write_tsv(.,file="silent_genes.bed",col_names = FALSE,quote = "none")
lowly_expressed_genes <- lowly_expressed[,c("chr","start","end","gene_id","PA_mean","strand","symbol")] %>% write_tsv(.,file="lowly_expressed_genes.bed",col_names = FALSE,quote = "none")
highly_expressed_genes <- highly_expressed[,c("chr","start","end","gene_id","PA_mean","strand","symbol")] %>% write_tsv(.,file="highly_expressed_genes.bed",col_names = FALSE,quote = "none")
```

# plot metagene agg plots for H3K27me1/2/3
```bash
# example
regions="silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed"
samples="PA_H3K27me1_merged.bw SETD2KO_H3K27me1_merged.bw K36MOE_H3K27me1_merged.bw DKO_H3K27me1_merged.bw TKO_H3K27me1_merged.bw"
sampleLabels="PA SETD2KO K36M-OE DKO TKO"
# cm
computeMatrix scale-regions -R ${regions} \
-S ${samples} \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
--missingDataAsZero --skipZeros \
-p 6 --samplesLabel ${sampleLabels} --verbose \
-o out.mat.gz
# plot
plotProfile -m out.mat.gz -o aggregate_profile.pdf --dpi 600 --colors blue hotpink saddlebrown gray darkkhaki --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --perGroup -y "MS-norm signal" --startLabel "TSS" --endLabel "TES" --legendLocation "best" --labelRotation 90 --plotFileFormat pdf
```

# plotting gene expression changes of genes with significant H3K27me3 spreading
```r
'%!in%' <- function(x,y)!('%in%'(x,y))
load("resLFC_TKO_PA.RData")
load("mMSC_RNASeq_dds.RData")
dds_RNAseq <- dds
rm(dds)
# log2FoldChanges 
dds_RNAseq$condition <- relevel(dds_RNAseq$condition,"NSD_unedit")
# model.matrix(~ condition
dds_RNAseq <- nbinomWaldTest(dds_RNAseq)
resultsNames_RNAseq <- resultsNames(dds_RNAseq) %>% as.data.frame()
# run shrinkage
resLFC_RNAseq <- lfcShrink(dds=dds_RNAseq,coef=4,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("id") %>% na.omit() %>% left_join(.,gene_symbol,by="id")
# filter by FPKM
counts_with_length <- fread("mm10_geneLengths.tsv") 
colnames(counts_with_length) <- c("id","length")
mcols(dds_RNAseq)$basepairs <- counts_with_length$length
norm.counts <- DESeq2::fpkm(dds_RNAseq,robust = TRUE) %>% as.data.frame() %>% tibble::rownames_to_column("id") %>% mutate(Mean_exp=(rowMeans(.[grep("TKO_rep1|TKO_rep2|PA_1|PA_2", names(.), value = TRUE)])),PA_mean_exp=(rowMeans(.[grep("PA_1|PA_2", names(.), value = TRUE)]))) %>% dplyr::select(c("id","Mean_exp","PA_mean_exp")) %>% dplyr::filter(Mean_exp > 1)
resLFC_RNAseq <- resLFC_RNAseq %>% mutate(quantile = ntile(baseMean,3)) %>% .[resLFC_RNAseq$id %in% norm.counts$id,]

```bash
featureCounts -a mm10_protein_coding_genes.saf -F SAF -o H3K27me3_mm10_protein_coding_genes.counts -T 6 -p *.bam
```

```r
raw_counts <- fread("H3K27me3_mm10_protein_coding_genes.counts") %>%
  as.data.frame() 

gene_length <- raw_counts$Length

mat <- raw_counts %>%
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()

metadata <- data.frame(kind = colnames(mat))

metadata$condition <- metadata$kind
# run DESeq2
metadata <- tibble::column_to_rownames(metadata, "kind")
dds <- DESeqDataSetFromMatrix(
  countData = mat,
  colData = metadata,
  design = ~condition
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

dds <- DESeq(dds)

resLFC_K27me3 <- lfcShrink(dds=dds,coef=2,type = "apeglm") %>% as.data.frame() %>% rownames_to_column("gene_name") %>% na.omit()

genes_with_increased_K27me3 <- resLFC_K27me3 %>% dplyr::filter(log2FoldChange > 0 & padj < 0.05)
# random genes
random_genes <- resLFC_K27me3 %>% .[.$gene_name %!in% genes_with_increased_K27me3$gene_name,]
RNAseq_for_genes_that_gained_K27me3 <- resLFC_RNAseq %>% .[.$gene_name %in% genes_with_increased_K27me3$gene_name,] %>% dplyr::select(log2FoldChange)
RNAseq_for_genes_that_gained_K27me3$cond <- "Genes with significant\nH3K27me3 invasion\nin their gene body"
nrow(RNAseq_for_genes_that_gained_K27me3 %>% dplyr::filter(log2FoldChange < 0))
set.seed(1)
RNAseq_for_random_genes <- resLFC_RNAseq %>% .[.$gene_name %in% random_genes$gene_name,] %>% dplyr::select(log2FoldChange) %>% dplyr::slice_sample(n=nrow(RNAseq_for_genes_that_gained_K27me3),replace=FALSE)
RNAseq_for_random_genes$cond <- "Random genes"
# combine
combined_FC <- rbind(RNAseq_for_genes_that_gained_K27me3,RNAseq_for_random_genes)
combined_FC$cond <- factor(combined_FC$cond,levels = c("Random genes","Genes with significant\nH3K27me3 invasion\nin their gene body"))
yaxis_label <- "Gene expression log2FoldChange (TKO/PA)"

second_group <- RNAseq_for_genes_that_gained_K27me3
first_group <- RNAseq_for_random_genes

# run stats
stats <- compare_means(log2FoldChange ~ cond,data=combined_FC,method="wilcox.test") %>% dplyr::mutate(pvalue=paste0("Wilcoxon, p=",.$p.format))
# run ggplot
ggplot(data = combined_FC,aes(x=cond,y=log2FoldChange,fill=cond)) +
  geom_violinhalf(show.legend = FALSE,position=position_nudge(x=0.1),size=0.1,linewidth=0.1) +
  geom_jitter(show.legend=FALSE,position=position_jitter(width=0.05, height=0.05),size=0.1,alpha=0.6) +
  geom_boxplot(width=0.2, color="black", alpha=0.5,outlier.size = 0.01,show.legend = FALSE,aes(colour="blue"),position=position_nudge(x=-0.2),linewidth=0.1) +
  scale_fill_manual(values=c("slateblue","lightcoral")) +
  geom_hline(yintercept =0,linetype="dashed",size=0.5,alpha=0.3,color="black",linewidth=0.3) +
  labs(x="",y="Gene expression log2 (FC) TKO/PA",title="PA versus TKO") +
  # stats
  geom_signif(textsize = 2,tip_length = 0.02,xmin=stats$group1,xmax=stats$group2,annotations=stats$p.signif,y_position = max(c(combined_FC$log2FoldChange * 1.05,combined_FC$log2FoldChange * 1.05,combined_FC$log2FoldChange*1.05)),size=0.5) +
  annotate("text",y=2, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange > 0,])),size=2) +
  annotate("text",y=2, x = 2.5, label = paste0(nrow(second_group[second_group$log2FoldChange > 0,])),size=2) +
  annotate("text",y=-2, x = 1.5, label = paste0(nrow(first_group[first_group$log2FoldChange < 0,])),size=2) +
  annotate("text",y=-2, x = 2.5, label = paste0(nrow(second_group[second_group$log2FoldChange < 0,])),size=2) +
  coord_flip() +
  theme(
    # panel
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y.left  = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.line.x.bottom  = element_line(colour="black"),
    axis.line.x.top  = element_line(colour="black"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed",linewidth=0.3),
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
ggsave(filename = "TKO_PA.K27me3_invasion.log2FoldChange_gene_expression.pdf",device = "pdf",units = "cm",width = 6.7,height=5,dpi = 600,bg="white")
```

# peakiness scores
```r
# sample script
bl <- import.bed("mm10_blacklist.bed")
rx <- fread("H3K27me3_rx.csv")
b <- fread("10T_PA_rep1_H3K27me3.1kb.bed",col.names = c('chr', 'start', 'end', 'score'))
samp <- "PA_1"
librarySize <- rx$chip[rx$samp == samp]
gr <- makeGRangesFromDataFrame(b,keep.extra.columns = TRUE)
gr$score <- (gr$score/librarySize)*1000000
top1 <- gr[gr$score > (quantile(gr$score,c(0.99)))]
top1_filt <- top1[!overlapsAny(top1,bl)]
p <- data.frame(samp=samp,
                score=mean(top1_filt$score))
# for subsequent samples
b <- fread("10T_PA_rep2_H3K27me3.1kb.bed",col.names = c('chr', 'start', 'end', 'score'))
samp <- "PA_2"
librarySize <- rx$chip[rx$samp == samp]
gr <- makeGRangesFromDataFrame(b,keep.extra.columns = TRUE)
gr$score <- (gr$score/librarySize)*1000000
top1 <- gr[gr$score > (quantile(gr$score,c(0.99)))]
top1_filt <- top1[!overlapsAny(top1,bl)]
s <- data.frame(samp=samp,
                score=mean(top1_filt$score))
p <- p %>% rows_insert(s)

# after running the previous scripts for all samples, then refactor the dataframe based on desired order of samples
p$condition <- factor(p$condition,levels=c("QuiKO","QKO","TKO","H3K36M-OE","DKO","SETD2KO","PA"))

stats <- aggregate(score ~ condition, p, function(x) c(mean = mean(x), sd = sd(x)))

s <- stats$score %>% as.data.frame() %>% mutate(condition = stats$condition)

ggplot() + stat_summary(mapping=aes(x=condition,y=score,fill=condition),data = p,geom="col",fun=mean,show.legend = FALSE) +
  geom_jitter(mapping=aes(x=condition,y=score,fill=condition),data = p,show.legend = FALSE,size=0.7) +
  geom_errorbar(mapping = aes(y=mean,x=condition,ymin = mean - sd, ymax = mean + sd),data = s, width = 0.2) +
  coord_flip() +
  scale_fill_manual(values=c("maroon3","orange","khaki","chocolate4","grey33","hotpink","blue")) +
  labs(x="",y="") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.y = element_text(family="Helvetica",size=9,colour = "black"),
    axis.text.x = element_text(family="Helvetica",size=9,colour = "black"),
    axis.text.y= element_text(family="Helvetica",size=9,colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    strip.text.x = element_text(size = 9),
    plot.margin=unit(c(0,0,0,0), "cm"))

ggsave(filename = "H3K27me3_peakiness_scores.pdf",path="figs",device = "pdf",dpi = 600,bg="white")
```

# agg plot centered on genes
```bash
computeMatrix scale-regions -R silent_genes.bed lowly_expressed_genes.bed highly_expressed_genes.bed \
-S mMSC_PA_merged_H3K27me3_msNorm.bw mMSC_QKO_merged_H3K27me3_msNorm.bw mMSC_QuiKO_merged_H3K27me3_msNorm.bw \
-bl mm10_blacklist.bed \
-b 2000 -a 2000 -m 20000 -bs 500 \
--missingDataAsZero --skipZeros \
-p 16 --samplesLabel PA QKO QuiKO --verbose \
-o ${out_mat}

plotProfile -m ${out_mat} -o ${out_png}.aggregate_profile.png --dpi 600 --colors blue orange maroon --numPlotsPerRow 3 --regionsLabel "silent" "lowly expressed" "highly expressed" --plotTitle "" --plotHeight 6.5 --plotWidth 8 --perGroup -y "MS-norm CPM" --startLabel "TSS" --endLabel "TES" --legendLocation "upper-right" --labelRotation 90
```